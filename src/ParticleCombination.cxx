#include "ParticleCombination.h"

#include "container_utils.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "logging.h"
#include "Model.h"

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
    daughter->Parent_ = shared_from_this();

    // add daughter to vector
    Daughters_.push_back(daughter);

    // copy daughter's indices into Indices_
    Indices_.insert(Indices_.end(), Daughters_.back()->indices().begin(), Daughters_.back()->indices().end());
}

//-------------------------
const std::shared_ptr<const ParticleCombination> ParticleCombination::origin() const
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
bool ParticleCombination::contains(const std::shared_ptr<const ParticleCombination>& B) const
{
    std::set<unsigned> setA(Indices_.begin(), Indices_.end());

    const std::vector<unsigned>& indicesB(B->indices());
    std::set<unsigned> setB(indicesB.begin(), indicesB.end());

    // check if B is subset of A
    return std::includes(setA.begin(), setA.end(), setB.begin(), setB.end());
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
bool disjoint(const ParticleCombinationVector& pcv)
{
    // check for any entry in pcv overlapping with another
    for (size_t i = 0; i < pcv.size(); ++i)
        for (size_t j = i + 1; j < pcv.size(); ++j)
            if (overlap(pcv[i]->indices(), pcv[j]->indices()))
                return false;
    return true;
}

//-------------------------
bool is_initial_state_particle_combination(const ParticleCombination& pc, const Model* m)
{
    auto p = pc.origin();

    for (auto& isp : m->initialStateParticles())
        if (any_of(isp.first->particleCombinations(), p))
            return true;

    return false;
}

//-------------------------
void prune_particle_combinations(ParticleCombinationVector& PCs)
{
    // get ISP PCs
    size_t ispNIndices(0);
    for (auto& pc : PCs)
        ispNIndices = std::max(ispNIndices, pc->origin()->indices().size());

    PCs.erase(std::remove_if(PCs.begin(), PCs.end(),
    [&] (std::shared_ptr<yap::ParticleCombination>& pc) { return (pc->origin()->indices().size() < ispNIndices); }),
    PCs.end());

    if (PCs.empty())
        throw exceptions::Exception("ParticleCombinations empty after pruning", "prune_particle_combinations");
}

//-------------------------
std::string indices_string(const ParticleCombination& pc, std::string before, std::string after)
{
    if (pc.indices().empty())
        return "(empty)";
    std::string s = before;
    for (auto i : pc.indices())
        s += std::to_string(i);
    s += after;
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
std::string to_string_with_parent(const ParticleCombination& pc)
{
    if (!pc.parent())
        return to_string(pc);

    return to_string(pc) + " in " + to_string(*pc.origin());
}

//-------------------------
// Comparison stuff:

ParticleCombination::Equal ParticleCombination::equalBySharedPointer;
ParticleCombination::EqualDown ParticleCombination::equalDown;
ParticleCombination::EqualUp ParticleCombination::equalUp;
ParticleCombination::EqualUpAndDown ParticleCombination::equalUpAndDown;
ParticleCombination::EqualByOrderedContent ParticleCombination::equalByOrderedContent;
ParticleCombination::EqualByOrderlessContent ParticleCombination::equalByOrderlessContent;
ParticleCombination::EqualDownByOrderlessContent ParticleCombination::equalDownByOrderlessContent;
ParticleCombination::EqualZemach ParticleCombination::equalZemach;

//-------------------------
bool ParticleCombination::EqualByOrderedContent::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
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
bool ParticleCombination::EqualDown::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
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
        if (!equalDown(A->daughters()[i], B->daughters()[i]))
            return false;

    // a match!
    return true;
}

//-------------------------
bool ParticleCombination::EqualUp::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equalByOrderedContent(A, B))
        return false;

    // if no more parents, return true
    if (!A->parent() and !B->parent())
        return true;

    // else continue up
    return ParticleCombination::equalUp(A->parent(), B->parent());
}

//-------------------------
bool ParticleCombination::EqualUpAndDown::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    if (!ParticleCombination::equalDown(A, B))
        return false;

    // if parents, return true
    if (!A->parent() and !B->parent())
        return true;

    // else check up
    return ParticleCombination::equalUp(A->parent(), B->parent());
}

//-------------------------
bool ParticleCombination::EqualByOrderlessContent::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
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
bool ParticleCombination::EqualDownByOrderlessContent::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!ParticleCombination::equalByOrderlessContent(A, B))
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size()) {
        return false;
    }

    unsigned matches(0);
    for (unsigned i = 0; i < A->daughters().size(); ++i)
        for (unsigned j = 0; j < B->daughters().size(); ++j)
            if (ParticleCombination::equalByOrderlessContent(A->daughters()[i], B->daughters()[j])) {
                ++matches;
                break;
            }

    if (matches == A->daughters().size())
        // a match!
        return true;

    return false;
}

//-------------------------
bool ParticleCombination::EqualByReferenceFrame::operator()(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B) const
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    // if both are nullptr, also return true
    if (A == B)
        return true;

    if (!ParticleCombination::equalByOrderlessContent(A->parent(), B->parent()))
        return false;

    return operator()(A->parent(), B->parent());
}

//-------------------------
bool ParticleCombination::EqualZemach::operator()(const std::shared_ptr<const ParticleCombination>& A,
        const std::shared_ptr<const ParticleCombination>& B) const
{
    //check if either empty
    if (!A or !B)
        return false;

    if (A->indices().size() > 3 or B->indices().size() > 3)
        throw exceptions::Exception("Zemach formalism cannot be used with 4 or more particles",
                                    "ParticleCombination::EqualZemach::operator()");

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
                                    "ParticleCombination::EqualZemach::operator()");
    auto rB = B->daughters()[0];
    auto sB = B->daughters()[1];
    if (sB->indices().size() == 2)
        std::swap(sB, rB);
    if (rB->indices().size() != 2 and sB->indices().size() != 1)
        throw exceptions::Exception("could not find resonance and spectator in B",
                                    "ParticleCombination::EqualZemach::operator()");

    return ParticleCombination::equalByOrderlessContent(rA, rB)
           and ParticleCombination::equalByOrderlessContent(sA, sB);
}

}
