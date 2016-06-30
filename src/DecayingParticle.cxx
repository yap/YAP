#include "DecayingParticle.h"

#include "BlattWeisskopf.h"
#include "CalculationStatus.h"
#include "Constants.h"
#include "container_utils.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "FreeAmplitude.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "SpinAmplitude.h"
#include "StatusManager.h"

#include <functional>
#include <iomanip>
#include <memory>

namespace yap {

//-------------------------
DecayingParticle::DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    Particle(q, mass, name),
    RadialSize_(std::make_shared<RealParameter>(radialSize))
{
}

//-------------------------
bool DecayingParticle::consistent() const
{
    bool C = Particle::consistent();

    if (RadialSize_->value() <= 0.) {
        FLOG(ERROR) << "Radial size not positive.";
        C &= false;
    }

    if (Channels_.empty()) {
        FLOG(ERROR) << "no channels specified.";
        return false; // further checks require at least one channel
    }

    // check no channel is empty
    if (std::any_of(Channels_.begin(), Channels_.end(), std::logical_not<DecayChannelVector::value_type>())) {
        FLOG(ERROR) << "DecayChannel vector contains nullptr";
        C &= false;
    }
    // check all channels' parents point to this
    if (std::any_of(Channels_.begin(), Channels_.end(), [&](const std::shared_ptr<DecayChannel>& dc) {return dc and dc->decayingParticle() != this;})) {
        FLOG(ERROR) << "DecayChannel vector contains channel not pointing back to this";
        C &= false;
    }
    // check consistency of all channels
    std::for_each(Channels_.begin(), Channels_.end(), [&](const std::shared_ptr<DecayChannel>& dc) {if (dc) C &= dc->consistent();});

    // check if all channels lead to same final state particles
    /// \todo This isn't necessary, we should think how to change this. Example: D -> KKpipi, with f0->KK and f0->pipi
    std::vector<std::shared_ptr<FinalStateParticle> > fsps0 = Channels_[0]->finalStateParticles();
    std::sort(fsps0.begin(), fsps0.end());
    for (unsigned i = 1; i < Channels_.size(); ++i) {
        std::vector<std::shared_ptr<FinalStateParticle> > fsps = Channels_[i]->finalStateParticles();
        std::sort(fsps.begin(), fsps.end());
        if (fsps != fsps0) {
            FLOG(ERROR) << "final state of channel " << i << " does not match.";
            C &= false;
        }
    }

    return C;
}

//-------------------------
std::shared_ptr<DecayChannel> DecayingParticle::addChannel(std::shared_ptr<DecayChannel> c)
{
    if (!c)
        throw exceptions::Exception("DecayChannel empty", "DecayingParticle::addChannel");

    if (c->particleCombinations().empty())
        throw exceptions::Exception(std::string("DecayChannel has no ParticleCombinations - ") + to_string(*c),
                                    "DecayingParticle::addChannel");

    // check ISP
    if (!Channels_.empty() and c->model() != model())
        throw exceptions::Exception("Model mismatch", "DecayingParticle::addChannel");

    // check if valid for DecayingParticle
    checkDecayChannel(c);

    Channels_.push_back(c);
    Channels_.back()->setDecayingParticle(this);

    // now that Model is set, register with Model (repeated registration has no effect)
    registerWithModel();

    // add particle combinations
    for (auto pc : Channels_.back()->particleCombinations())
        addParticleCombination(pc);

    /// create decay trees for channel:

    /// loop over spin amplitudes of channel
    for (auto& sa : Channels_.back()->spinAmplitudes()) {

        // loop over possible parent spin projections
        for (const auto& two_M : sa->twoM()) {

            // create DecayTree with new FreeAmplitude
            auto DT_M = DecayTree(std::make_shared<FreeAmplitude>(Channels_.back(), sa, two_M));
            // add BlattWeisskopf object to it
            modifyDecayTree(DT_M);

            // loop over possible daughter spin projections
            for (const auto& two_m : sa->twoM(two_M)) {

                // loop over particle combinations of this decaying particle
                for (const auto& pc : particleCombinations()) {

                    // initialize vector of possible decay trees from
                    // initial DecayTree created above
                    DecayTreeVector DTV(1, std::make_shared<DecayTree>(DT_M));

                    // loop over daughters in channel
                    for (size_t d = 0; d < Channels_.back()->daughters().size(); ++d) {

                        // try to cast daughter to decaying particle
                        auto dp = std::dynamic_pointer_cast<DecayingParticle>(Channels_.back()->daughters()[d]);

                        // if decaying particle
                        if (dp) {

                            // check if daughter has any decay trees with spin projection
                            if (dp->DecayTrees_.find(two_m[d]) != dp->DecayTrees_.end()) {

                                // create temp tree vector to store new copies into
                                DecayTreeVector DTV_temp;

                                // loop over decay trees of daughter with appropriate spin projection
                                for (const auto& dt : dp->DecayTrees_[two_m[d]]) {
                                    // check that decay channel of free amplitude of decay tree has particle combination
                                    if (any_of(dt->freeAmplitude()->decayChannel()->particleCombinations(), pc->daughters()[d])) {
                                        for (const auto& DT : DTV) {
                                            // add copy of DT to DTV_temp
                                            DTV_temp.push_back(std::make_shared<DecayTree>(*DT));
                                            // add decay tree to it
                                            DTV_temp.back()->setDaughterDecayTree(d, dt);
                                        }
                                    }
                                }

                                // replace DTV with DTV_temp
                                DTV = DTV_temp;
                            }
                            // else clear DTV
                            else
                                DTV.clear();

                        }
                        // else not decaying particle
                        else {

                            // check that particle has daughter particle combination
                            if (any_of(Channels_.back()->daughters()[d]->particleCombinations(), pc->daughters()[d])) {
                                // set all elts of DTV to have proper daughter spin projection
                                for (auto& DT : DTV)
                                    DT->setDaughterSpinProjection(d, two_m[d]);
                            }
                            // else clear DTV
                            else
                                DTV.clear();
                        }

                        // if DTV now empty, break
                        if (DTV.empty())
                            break;

                    } // ends loop over daughters

                    // if decay trees were created, add them into DecayTrees_
                    // if they aren't already present in it
                    if (!DTV.empty()) {
                        auto& dtv_M = DecayTrees_[two_M];
                        for (const auto& DT : DTV) {
                            // check that DT isn't already equal to one in dtv_M
                            if (std::none_of(dtv_M.begin(), dtv_M.end(), [&DT](const std::shared_ptr<DecayTree>& dt) {return *dt == *DT;}))
                            dtv_M.push_back(DT);
                        }
                    }

                } // ends loop over particle combinations of this decaying particle
            } // ends loop over spin projections of daughters
        } // ends loop over spin projection of parent
    } // ends loop over spin amplitude

    return Channels_.back();
}

//-------------------------
std::shared_ptr<DecayChannel> DecayingParticle::addChannel(const ParticleVector& daughters)
{
    return addChannel(std::make_shared<DecayChannel>(daughters));
}

//-------------------------
const Model* DecayingParticle::model() const
{
    return Channels_.empty() ? nullptr : Channels_[0]->model();
}

//-------------------------
void DecayingParticle::addParticleCombination(const std::shared_ptr<ParticleCombination>& pc)
{
    Particle::addParticleCombination(pc);

    // add also to all BlattWeiskopf barrier factors
    for (auto& kv : BlattWeisskopfs_)
        if (pc->daughters().size() == 2)
            kv.second->addParticleCombination(pc);

    // add to DecayChannels,
    // if DecayChannel contains particle combination with same content (without checking parent)
    // this is for the setting of ParticleCombination's with parents
    for (auto& dc : Channels_) {
        if (any_of(dc->particleCombinations(), pc, ParticleCombination::equalDown))
            dc->addParticleCombination(pc);
    }
}

//-------------------------
void DecayingParticle::pruneParticleCombinations()
{
    Particle::pruneParticleCombinations();

    for (auto& c : Channels_)
        c->pruneParticleCombinations();
}

//-------------------------
void DecayingParticle::fixSolitaryFreeAmplitudes()
{
    // loop over entries in map of (spin projection) -> (decay tree vector)
    for (auto& m_dtv : DecayTrees_)
        // if only available decay tree
        if (m_dtv.second.size() == 1)
            m_dtv.second[0]->freeAmplitude()->setVariableStatus(VariableStatus::fixed);
    for (auto& c : Channels_)
        c->fixSolitaryFreeAmplitudes();
}

//-------------------------
std::vector< std::shared_ptr<FinalStateParticle> > DecayingParticle::finalStateParticles(unsigned i) const
{
    if (i >= Channels_.size())
        throw exceptions::Exception("Channel index too high (" + std::to_string(i) + " >= " + std::to_string(Channels_.size()) + ")",
                                    "DecayingParticle::finalStateParticles");
    if (!Channels_[i])
        throw exceptions::Exception("Channel " + std::to_string(i) + " is empty", "DecayingParticle::finalStateParticles");

    return Channels_[i]->finalStateParticles();
}

//-------------------------
std::string DecayingParticle::decayChainLevelAsString(int level) const
{
    std::ostringstream os;

    // get maximum length of particle names
    static size_t padding = 0;
    static size_t paddingSpinAmp = 0;
    if (padding == 0 || level == -1) {
        padding = std::max(padding, name().length());
        for (auto& c : Channels_) {
            for (std::shared_ptr<Particle> d : c->daughters()) {
                padding = std::max(padding, d->name().length());
                if (std::dynamic_pointer_cast<DecayingParticle>(d))
                    os << std::static_pointer_cast<DecayingParticle>(d)->decayChainLevelAsString(-1);
                paddingSpinAmp = std::max(paddingSpinAmp, to_string(c->spinAmplitudes()).length());
            }
        }
        if (level == -1)
            return os.str();
    }

    for (size_t i = 0; i < Channels_.size(); ++i) {
        if (i > 0)
            os << "\n" << std::setw(level * (padding * 3 + 8 + paddingSpinAmp)) << "";

        os << std::left << std::setw(padding) << this->name() << " ->";
        for (std::shared_ptr<Particle> d : Channels_[i]->daughters())
            os << " " << std::setw(padding) << d->name();
        os << std::left << std::setw(paddingSpinAmp)
                  << to_string(Channels_[i]->spinAmplitudes());

        for (std::shared_ptr<Particle> d : Channels_[i]->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                os << ",  ";
                os << std::static_pointer_cast<DecayingParticle>(d)->decayChainLevelAsString(level + 1);
            }
    }

    if (level == 0)
        os << "\n";

    return os.str();
}

//-------------------------
std::shared_ptr<DecayChannel> DecayingParticle::channel(const ParticleVector& daughters)
{
    auto it = std::find_if(Channels_.begin(), Channels_.end(), [&](std::shared_ptr<DecayChannel> dc) {return dc->daughters() == daughters;});
    if (it == Channels_.end())
        throw exceptions::Exception("Channel not found", "DecayingParticle::channel");
    return *it;
}

//-------------------------
FreeAmplitudeSet DecayingParticle::freeAmplitudes() const
{
    FreeAmplitudeSet S;
    for (auto m_dtv : DecayTrees_)
        for (auto dt : m_dtv.second)
            if (dt->freeAmplitude())
                S.insert(dt->freeAmplitude());
    return S;
}

//-------------------------
FreeAmplitudeSet freeAmplitudes(const std::map<int, DecayTreeVector>& m_dtv_map)
{
    FreeAmplitudeSet S;
    for (auto& m_dtv : m_dtv_map) {
        auto s = freeAmplitudes(m_dtv.second);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
void DecayingParticle::storeBlattWeisskopf(unsigned l)
{
    // if BW is not already stored for L, add it
    if (BlattWeisskopfs_.find(l) == BlattWeisskopfs_.end())
        BlattWeisskopfs_.emplace(l, std::make_shared<BlattWeisskopf>(l, this));
}

//-------------------------
void DecayingParticle::modifyDecayTree(DecayTree& dt) const
{
    if (!dt.freeAmplitude())
        throw exceptions::Exception("DecayTree has nullptr free amplitude", "DecayingParticle::modifyDecayTree");

    if (!dt.freeAmplitude()->spinAmplitude())
        throw exceptions::Exception("FreeAmplitude's SpinAmplitude is nullptr", "DecayingParticle::modifyDecayTree");

    // find BlattWeisskopf object
    if (dt.freeAmplitude()->spinAmplitude()->L() > 0) {
        auto bw = BlattWeisskopfs_.find(dt.freeAmplitude()->spinAmplitude()->L());
        if (bw == BlattWeisskopfs_.end())
            throw exceptions::Exception("No Blatt-Weisskopf factor found for L = "
                                        + std::to_string(dt.freeAmplitude()->spinAmplitude()->L()),
                                        "DecayingParticle::modifyDecayTree");

        if (!bw->second)
            throw exceptions::Exception("BlattWeisskopf is nullptr", "DecayingParticle::modifyDecayTree");

        // Add BlattWeisskopf object
        dt.addDataAccessor(*bw->second);
    }
}

//-------------------------
const std::complex<double> amplitude(const std::map<int, DecayTreeVector>& m_dtv_map, const DataPoint& d)
{
    auto A = Complex_0;
    for (const auto& m_dtv : m_dtv_map)
        for (const auto& dt : m_dtv.second)
            A += amplitude(*dt, d);
    return A;
}


}
