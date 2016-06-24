#include "DecayTree.h"

#include "DecayChannel.h"
#include "Exceptions.h"
#include "FreeAmplitude.h"
#include "logging.h"
#include "RecalculableDataAccessor.h"
#include "SpinAmplitude.h"
#include "VariableStatus.h"

#include <iterator>
#include <functional>

namespace yap {

//-------------------------
DecayTree::DecayTree(std::shared_ptr<FreeAmplitude> free_amp) :
    FreeAmplitude_(free_amp),
    DaughtersTwoM_(FreeAmplitude_->spinAmplitude()->finalTwoJ().size()),
    FreeAmplitude_(free_amp)
{}

//-------------------------
const Model* DecayTree::model() const
{
    return (FreeAmplitude_) ? FreeAmplitude_->model() : nullptr;
}

//-------------------------
const std::complex<double> amplitude(const DecayTreeVector& dtv, const DataPoint& d)
{
    auto A = Complex_0;
    for (const auto& dt : dtv)
        A += amplitude(*dt, d);
    return A;
}

//-------------------------
FreeAmplitudeSet freeAmplitudes(const DecayTree& DT)
{
    FreeAmplitudeSet S = {DT.freeAmplitude()};
    for (auto& d_dt : DT.daughterDecayTrees()) {
        auto s = freeAmplitudes(*d_dt.second);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
FreeAmplitudeSet freeAmplitudes(const DecayTreeVector& DTV)
{
    FreeAmplitudeSet S;
    for (auto& DT : DTV) {
        auto s = freeAmplitudes(*DT);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
const std::complex<double> DecayTree::dataDependentAmplitude(const DataPoint& d) const
{
    auto A = Complex_0;
    FDEBUG("decayChannel " << to_string(*FreeAmplitude_->decayChannel()));
    for (const auto& pc : FreeAmplitude_->decayChannel()->particleCombinations()) {
        FDEBUG("pc " << to_string(*pc));
        A += dataDependentAmplitude(d, pc);
    }
    return A;
}

//-------------------------
const std::complex<double> DecayTree::dataDependentAmplitude(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    FDEBUG("decayChannel " << to_string(*FreeAmplitude_->decayChannel()));
    FDEBUG("pc" << to_string(*pc));

    // spin amplitude
    auto A = FreeAmplitude_->spinAmplitude()->amplitude(d, pc, FreeAmplitude_->twoM(), DaughtersTwoM_[0], DaughtersTwoM_[1]);

    // recalculable amplitude
    for (const auto& rda : RecalculableDataAccessors_)
        A *= rda->value(d, pc);

    // likewise for daughters
    for (const auto& d_dt : DaughterDecayTrees_)
        A *= d_dt.second->dataDependentAmplitude(d, pc->daughters()[d_dt.first]);

    return A;
}

//-------------------------
const VariableStatus DecayTree::dataDependentAmplitudeStatus() const
{
    // check status of recalculable components, terminating as soon as one is seen to have changed
    for (const auto& rda : RecalculableDataAccessors_)
        if (variableStatus(*rda) == VariableStatus::changed)
            return VariableStatus::changed;
    // check daughters, terminating as soon as one is seen to have changed
    for (const auto& d_dt : DaughterDecayTrees_)
        if (d_dt.second->dataDependentAmplitudeStatus() == VariableStatus::changed)
            return VariableStatus::changed;
    // else unchanged
    return VariableStatus::unchanged;
}

//-------------------------
const std::complex<double> DecayTree::dataIndependentAmplitude() const
{
    // spin amplitude
    auto A = FreeAmplitude_->value();
    // likewise for daughters
    for (const auto& d_dt : DaughterDecayTrees_)
        A *= d_dt.second->dataIndependentAmplitude();
    return A;
}

//-------------------------
void DecayTree::setDaughterDecayTree(unsigned i, std::shared_ptr<DecayTree> dt)
{
    if (i >= DaughtersTwoM_.size())
        throw exceptions::Exception("index exceeds number of daughters", "DecayTree::setDaughterDecayTree");

    if (!dt->freeAmplitude())
        throw exceptions::Exception("DecayTree's free amplitude is nullptr", "DecayTree:setDaughterDecayTree");

    if (DaughterDecayTrees_.find(i) != DaughterDecayTrees_.end())
        throw exceptions::Exception("DecayTree for this daughter already set", "DecayTree::setDaughterDecayTree");

    DaughterDecayTrees_[i] = dt;
    setDaughterSpinProjection(i, dt->freeAmplitude()->twoM());
}

//-------------------------
void DecayTree::addDataAccessor(const RecalculableDataAccessor& rda)
{
    if (!FreeAmplitude_)
        throw exceptions::Exception("FreeAmplitude is nullptr", "DecayTree::addDataAccessor");

    if (!FreeAmplitude_->checkParticleCombinations(rda))
        throw exceptions::Exception("RecalculableDataAccessor doesn't have all ParticleCombinations required by FreeAmplitude",
                                    "DecayTree::addDataAccessor");

    RecalculableDataAccessors_.push_back(&rda);
}

//-------------------------
std::string DecayTree::asString(std::string offset) const
{
    auto s = to_string(*FreeAmplitude_);
    offset += "    ";
    for (const auto& d_dt : DaughterDecayTrees_)
        s += "\n" + offset + "d[" + std::to_string(d_dt.first) + "] --> " + d_dt.second->asString(offset);
    return s;
}

//-------------------------
unsigned depth(const DecayTree& DT)
{
    unsigned d = 0;
    for (const auto& i_dt : DT.daughterDecayTrees())
        d = std::max(d, depth(*i_dt.second));
    return d + 1;
}

//-------------------------
bool has_changed(const DecayTreeVector::value_type& dt)
{ return dt->dataDependentAmplitudeStatus() == VariableStatus::changed; }

//-------------------------
const DecayTreeVector select_changed(const DecayTreeVector& dtv)
{
    DecayTreeVector C;
    C.reserve(dtv.size());
    std::copy_if(dtv.begin(), dtv.end(), std::back_inserter(C),
                 std::function<bool(const DecayTreeVector::value_type&)>(has_changed));
    return C;
}

}
