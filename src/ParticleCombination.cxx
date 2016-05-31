#include "ParticleCombination.h"

#include "container_utils.h"
#include "Exceptions.h"
#include "logging.h"
#include "Model.h"
#include "Particle.h"

#include <algorithm>
#include <set>

namespace yap {

//-------------------------
void ParticleCombination::addDaughter(std::shared_ptr<ParticleCombination> daughter)
{
    if (daughter->indices().empty())
        throw exceptions::Exception("daughter contains no indices", "ParticleCombination::addDaughter");

    /// Check that new daughter does not share content with other daughters?
    if (overlap(daughter->indices(), Indices_))
        throw exceptions::Exception("daughter overlaps with other daughters", "ParticleCombination::addDaughter");

    /// Check that new daughter doesn't already have parent set
    if (daughter->parent())
        throw exceptions::Exception("daughter's parent is already set", "ParticleCombination::addDaughter");

    // set daughter's parent to shared_from_this
    std::const_pointer_cast<ParticleCombination>(daughter)->Parent_ = shared_from_this();

    // add daughter to vector
    Daughters_.push_back(daughter);

    // copy daughter's indices into Indices_
    Indices_.insert(Indices_.end(), Daughters_.back()->indices().begin(), Daughters_.back()->indices().end());
}

//-------------------------
std::shared_ptr<ParticleCombination> ParticleCombination::origin()
{
    auto pc = shared_from_this();
    while (pc->parent())
        pc = pc->parent();
    return pc;
}

//-------------------------
ParticleCombinationVector ParticleCombination::leaves()
{
    if (Daughters_.empty())
        return ParticleCombinationVector(1, shared_from_this());

    ParticleCombinationVector V;
    for (const auto& d : Daughters_) {
        auto v = d->leaves();
        V.insert(V.end(), v.begin(), v.end());
    }
    return V;
}

//-------------------------
bool ParticleCombination::decaysToFinalStateParticles() const
{
    for (auto& leaf : const_cast<ParticleCombination*>(this)->leaves())
        if (!leaf->isFinalStateParticle())
            return false;
    return true;
}

//-------------------------
std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > combinations(std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > >& P, Model* m)
{
    if (P.size() == 1)
        return P;

    DEBUG("construct combinations of: " << to_string(P));

    std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > rest;
    rest.insert(rest.end(), P.begin() + 1, P.end());

    auto combiRest = combinations(rest, m);
    return combinations(P[0], combiRest, m);
}

//-------------------------
std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > combinations(std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> >& A, std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > B, Model* m)
{
    std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > > result;

    DEBUG("A: " << to_string(A));
    DEBUG("B: " << to_string(B));

    for (auto& PCA : A) {
        for (auto& PCBs : B) {

            for (size_t i = 0; i < PCBs.size(); ++i) {
                auto& PCB = PCBs[i];

                std::vector<unsigned> indicesA;
                for (auto& pc : PCA.second)
                    indicesA.insert(indicesA.end(), pc->indices().begin(), pc->indices().end());

                std::vector<unsigned> indicesB;
                for (auto& pc : PCB.second)
                    indicesB.insert(indicesB.end(), pc->indices().begin(), pc->indices().end());


                // check that PCA and PCB don't overlap in FSP content
                if (overlap(indicesA, indicesB)) {

                    // make sure not to symmetrize for identical particles
                    if (PCA.first and PCA.first == PCB.first) {
                        DEBUG("same - " << to_string(PCA.second) << " - " << to_string(PCB.second));

                        PCBs.erase(std::find(PCBs.begin(), PCBs.end(), PCA));

                        --i;
                    }

                    continue;
                }

                ParticleCombinationVector AB;
                for (auto& pc : PCA.second)
                    AB.push_back(pc);

                for (auto& pc : PCB.second)
                    AB.push_back(pc);

                std::shared_ptr<Particle> part = PCA.first;
                if (part != PCB.first)
                    part.reset();

                result.push_back({std::make_pair(part, AB)});
            }
        }
    }

    return result;
}

//-------------------------
bool ParticleCombination::consistent() const
{
    bool C = true;

    // should have no daughters or 2 or more daughters:
    if (Daughters_.size() == 1) {
        FLOG(ERROR) << "has only one daughter.";
        C &= false;
    }

    if (Indices_.empty()) {
        FLOG(ERROR) << "has no indices.";
        C &= false;
    }

    // if has daughters
    if (!Daughters_.empty()) {
        // Check Indices_ doesn't have duplicates
        // create unique_copy of Indices_ (as set)
        std::set<unsigned> U(Indices_.begin(), Indices_.end());
        // check unique_copy-object's size == this object's size
        if (U.size() != Indices_.size()) {
            FLOG(ERROR) << "index vector contains duplicate entries (" << U.size() << " != " << Indices_.size() << ").";
            C &= false;
        }

        // count number of daughters with parent not set to this
        auto n = std::count_if(Daughters_.begin(), Daughters_.end(),
        [&](const ParticleCombinationVector::value_type & d) {return d->parent().get() != this;});
        if (n != 0) {
            FLOG(ERROR) << n << " daughters' parent not set to this ParticleCombination.";
            C &= false;
        }

        // check consistency of daughters
        std::for_each(Daughters_.begin(), Daughters_.end(),
        [&](const ParticleCombinationVector::value_type & d) {C &= d->consistent();});
    }
    // if Daugthers_ empty, should have one and only index (as FSP)
    else if (Indices_.size() != 1) {
        FLOG(ERROR) << "contains wrong number of indices for final-state particle (" << Indices_.size() << " != 1)";
        C &= false;
    }

    return C;
}

//-------------------------
std::string indices_string(const ParticleCombination& pc)
{
    if (pc.indices().empty())
        return "(empty)";
    std::string s = "(";
    for (auto i : pc.indices())
        s += std::to_string(i);
    s += ")";
    return s;
}

//-------------------------
std::string to_string(const ParticleCombination& pc)
{
    auto s = indices_string(pc);

    if (pc.daughters().empty())
        return s;

    s += " -> ";

    for (auto& d : pc.daughters()) {
        s += "(";
        for (auto i : d->indices())
            s += std::to_string(i);
        s += ") + ";
    }
    s.erase(s.size() - 3, 3);
    for (auto& d : pc.daughters())
        if (!d->isFinalStateParticle())
            s += "; " + to_string(*d);
    return s;
}


//-------------------------
std::string to_string(const ParticleCombinationVector& PCs)
{
    std::string s = "[";
    for (auto& pc : PCs) {
        s += to_string(*pc);
        s += "; ";
    }

    s.erase(s.size() - 2, 2);
    s += "]";

    return s;
}

//-------------------------
std::string to_string(const std::vector<ParticleCombinationVector>& PCs)
{
    std::string s = "{";
    for (auto& pc : PCs) {
        s += to_string(pc);
        s += "; ";
    }

    s.erase(s.size() - 2, 2);
    s += "}";

    return s;
}

//-------------------------
std::string to_string(const std::vector<std::vector<ParticleCombinationVector>>& PCs)
{
    std::string s = "{ ";
    for (auto& pc : PCs) {
        s += to_string(pc);
        s += "\n  ";
    }

    s.erase(s.size() - 3, 3);
    s += " }";

    return s;
}

//-------------------------
std::string to_string(const std::pair<std::shared_ptr<Particle>, ParticleCombinationVector>& pc)
{
    std::string s = (pc.first) ? pc.first->name() : "nullptr";
    s += ": ";
    s += to_string(pc.second);

    return s;
}

//-------------------------
std::string to_string(const std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> >& PCs)
{
    std::string s = "{";
    for (auto& pc : PCs) {
        s += to_string(pc);
        s += "; ";
    }

    s.erase(s.size() - 2, 2);
    s += "}";

    return s;
}

//-------------------------
std::string to_string(const std::vector<std::vector<std::pair<std::shared_ptr<Particle>, ParticleCombinationVector> > >& PCs)
{
    std::string s = "{ ";
    for (auto& pc : PCs) {
        s += to_string(pc);
        s += "\n  ";
    }

    s.erase(s.size() - 3, 3);
    s += " }";

    return s;
}

//-------------------------
// Comparison stuff:

ParticleCombination::Equiv ParticleCombination::equivBySharedPointer;
ParticleCombination::EquivDown ParticleCombination::equivDown;
ParticleCombination::EquivUp ParticleCombination::equivUp;
ParticleCombination::EquivUpAndDown ParticleCombination::equivUpAndDown;
ParticleCombination::EquivByOrderedContent ParticleCombination::equivByOrderedContent;
ParticleCombination::EquivByOrderlessContent ParticleCombination::equivByOrderlessContent;
ParticleCombination::EquivDownByOrderlessContent ParticleCombination::equivDownByOrderlessContent;
ParticleCombination::EquivZemach ParticleCombination::equivZemach;

//-------------------------
bool ParticleCombination::EquivByOrderedContent::operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // Check indices
    if (A->indices() != B->indices())
        return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivDown::operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // check ordered content
    if (A->indices() != B->indices())
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size())
        return false;

    for (unsigned i = 0; i < A->daughters().size(); ++i)
        if (!equivDown(A->daughters()[i], B->daughters()[i]))
            return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EquivUp::operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivByOrderedContent(A, B))
        return false;

    // if no more parents, return true
    if (!A->parent() and !B->parent())
        return true;

    // else continue up
    return ParticleCombination::equivUp(A->parent(), B->parent());
}

//-------------------------
bool ParticleCombination::EquivUpAndDown::operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const
{
    if (!ParticleCombination::equivDown(A, B))
        return false;

    // if parents, return true
    if (!A->parent() and !B->parent())
        return true;

    // else check up
    return ParticleCombination::equivUp(A->parent(), B->parent());
}

//-------------------------
bool ParticleCombination::EquivByOrderlessContent::operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // check size of indices vectors
    if (A->indices().size() != B->indices().size())
        return false;

    // check contents of indices vectors
    // (creating a set will sort entries for easy comparison,
    // since order doesn't matter)
    std::set<unsigned> a(A->indices().begin(), A->indices().end());
    std::set<unsigned> b(B->indices().begin(), B->indices().end());

    return std::equal(a.begin(), a.end(), b.begin());
}

//-------------------------
bool ParticleCombination::EquivDownByOrderlessContent::operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equivByOrderlessContent(A, B))
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size()) {
        return false;
    }

    unsigned matches(0);
    for (unsigned i = 0; i < A->daughters().size(); ++i)
        for (unsigned j = 0; j < B->daughters().size(); ++j)
            if (ParticleCombination::equivByOrderlessContent(A->daughters()[i], B->daughters()[j])) {
                ++matches;
                break;
            }

    if (matches == A->daughters().size())
        // a match!
        return true;

    return false;
}

//-------------------------
bool ParticleCombination::EquivByReferenceFrame::operator()(const std::shared_ptr<ParticleCombination>& A, const std::shared_ptr<ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    // if both are nullptr, also return true
    if (A == B)
        return true;

    if (!ParticleCombination::equivByOrderlessContent(A->parent(), B->parent()))
        return false;

    return operator()(A->parent(), B->parent());
}

//-------------------------
bool ParticleCombination::EquivZemach::operator()(const std::shared_ptr<ParticleCombination>& A,
        const std::shared_ptr<ParticleCombination>& B) const
{
    //check if either empty
    if (!A or !B)
        return false;

    if (A->indices().size() > 3 or B->indices().size() > 3)
        throw exceptions::Exception("Zemach formalism cannot be used with 4 or more particles",
                                    "ParticleCombination::EquivZemach::operator()");

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // check if both size < 3
    if (A->indices().size() < 3 and B->indices().size() < 3)
        return true;

    // now check if sizes same
    if (A->indices().size() != B->indices().size())
        return false;

    // find resonance and spectator
    auto rA = A->daughters()[0];
    auto sA = A->daughters()[1];
    if (sA->indices().size() == 2)
        std::swap(sA, rA);
    if (rA->indices().size() != 2 and sA->indices().size() != 1)
        throw exceptions::Exception("could not find resonance and spectator in A",
                                    "ParticleCombination::EquivZemach::operator()");
    auto rB = B->daughters()[0];
    auto sB = B->daughters()[1];
    if (sB->indices().size() == 2)
        std::swap(sB, rB);
    if (rB->indices().size() != 2 and sB->indices().size() != 1)
        throw exceptions::Exception("could not find resonance and spectator in B",
                                    "ParticleCombination::EquivZemach::operator()");

    return ParticleCombination::equivByOrderlessContent(rA, rB)
           and ParticleCombination::equivByOrderlessContent(sA, sB);
}

}
