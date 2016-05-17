#include "DecayingParticle.h"

#include "AmplitudePair.h"
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
std::complex<double> DecayingParticle::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, int two_m, StatusManager& sm) const
{
    auto A = Complex_0;

    // sum up DecayChannel::amplitude over each channel
    for (const auto& c : channels())
        if (c->hasParticleCombination(pc))
            A += c->amplitude(d, pc, two_m, sm);

    return A;
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
    if (std::any_of(Channels_.begin(), Channels_.end(), [](const std::shared_ptr<DecayChannel>& dc) {return !dc;})) {
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
    std::vector<std::shared_ptr<FinalStateParticle> > fsps0 = finalStateParticles(0);
    std::sort(fsps0.begin(), fsps0.end());
    for (unsigned i = 1; i < nChannels(); ++i) {
        std::vector<std::shared_ptr<FinalStateParticle> > fsps = finalStateParticles(i);
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

    Channels_.push_back(c);
    Channels_.back()->setDecayingParticle(this);

    // now that Model is set, register with Model (repeated registration has no effect)
    registerWithModel();

    // if this is to be the initial state particle
    if (!model()->initialStateParticle() and finalStateParticles().size() == model()->finalStateParticles().size())
        const_cast<Model*>(static_cast<const DecayingParticle*>(this)->model())->setInitialStateParticle(std::static_pointer_cast<DecayingParticle>(shared_from_this()));

    // add particle combinations
    for (auto pc : Channels_.back()->particleCombinations())
        addParticleCombination(pc);

    /// create decay trees for channel:

    /// loop over spin amplitudes of channel
    for (auto& sa_apm : Channels_.back()->Amplitudes_) {

        // loop over possible parent spin projections
        for (const auto& M_m : sa_apm.first->amplitudes()) {

            // create DecayTree with new FreeAmplitude
            auto DT_M = DecayTree(std::make_shared<FreeAmplitude>(Channels_.back(), sa_apm.first, M_m.first));
            // add BlattWeisskopf object to it
            modifyDecayTree(DT_M);

            // loop over possible daughter spin projections
            for (const auto& m_cdv : M_m.second) {

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

                            // create temp tree vector to store new copies into
                            DecayTreeVector DTV_temp;

                            // check if daughter has any decay trees with spin projection
                            if (dp->DecayTrees_.find(m_cdv.first[d]) != dp->DecayTrees_.end()) {
                                // loop over decay trees of daughter with appropriate spin projection
                                for (const auto& dt : dp->DecayTrees_[m_cdv.first[d]]) {
                                    // check that decay channel of free amplitude of decay tree has particle combination
                                    if (dt->freeAmplitude()->decayChannel()->hasParticleCombination(pc->daughters()[d])) {
                                        for (const auto& DT : DTV) {
                                            // add copy of DT to DTV_temp
                                            DTV_temp.push_back(std::make_shared<DecayTree>(*DT));
                                            // add decay tree to it
                                            DTV_temp.back()->setDaughterDecayTree(d, dt);
                                        }
                                    }
                                }
                            }
                            // else clear DTV
                            else {
                                DTV.clear();
                            }

                            // replace DTV with DTV_temp
                            DTV = DTV_temp;
                        }
                        // else non-decaying particle
                        else {

                            // check that particle has daughter particle combination
                            if (Channels_.back()->daughters()[d]->hasParticleCombination(pc->daughters()[d])) {
                                // set all elts of DTV to have proper daughter spin projection
                                for (auto& DT : DTV)
                                    DT->setDaughterSpinProjection(d, m_cdv.first[d]);
                            }
                            // else clear DTV
                            else {
                                DTV.clear();
                            }
                        }

                        // if DTV now empty, break
                        if (DTV.empty())
                            break;

                    } // ends loop over daughters

                    // if decay trees were created, add them into DecayTrees_
                    // if they aren't already present in it
                    if (!DTV.empty()) {
                        auto& dtv_M = DecayTrees_[M_m.first];
                        for (const auto& DT : DTV) {
                            // check that DT isn't already equivalent to one in dtv_M
                            if (std::none_of(dtv_M.begin(), dtv_M.end(), [&DT](const std::shared_ptr<DecayTree>& dt) {return *dt == *DT;}))
                            dtv_M.push_back(DT);
                        }
                    }

                } // ends loop over particle combinations of this decaying particle
            } // ends loop over spin projections of daughters
        } // ends loop over spin projection of parent
    } // ends loop over spin amplitude

    FDEBUG(*Channels_.back() << " with N(PC) = " << Channels_.back()->particleCombinations().size());
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
void DecayingParticle::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    Particle::addParticleCombination(pc);

    // add also to all BlattWeiskopf barrier factors
    for (auto& kv : BlattWeisskopfs_)
        kv.second->addParticleCombination(pc);

    // add to DecayChannels,
    // if DecayChannel contains particle combination with same content (without checking parent)
    // this is for the setting of ParticleCombination's with parents
    for (auto& dc : Channels_) {
        if (dc->hasParticleCombination(pc, ParticleCombination::equivDown))
            dc->addParticleCombination(pc);
    }

    // check if also model's initial state particle
    if (model() and model()->initialStateParticle() == shared_from_this())
        const_cast<Model*>(static_cast<const DecayingParticle*>(this)->model())->addParticleCombination(pc);
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
void DecayingParticle::printDecayChainLevel(int level) const
{
    // get maximum length of particle names
    static size_t padding = 0;
    static size_t paddingSpinAmp = 0;
    if (padding == 0 || level == -1) {
        padding = std::max(padding, name().length());
        for (auto& c : Channels_) {
            for (std::shared_ptr<Particle> d : c->daughters()) {
                padding = std::max(padding, d->name().length());
                if (std::dynamic_pointer_cast<DecayingParticle>(d))
                    std::static_pointer_cast<DecayingParticle>(d)->printDecayChainLevel(-1);
                paddingSpinAmp = std::max(paddingSpinAmp, to_string(c->spinAmplitudes()).length());
            }
        }
        if (level == -1)
            return;
    }

    for (unsigned int i = 0; i < nChannels(); ++i) {
        if (i > 0)
            std::cout << "\n" << std::setw(level * (padding * 3 + 8 + paddingSpinAmp)) << "";

        std::cout << std::left << std::setw(padding) << this->name() << " ->";
        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            std::cout << " " << std::setw(padding) << d->name();
        std::cout << std::left << std::setw(paddingSpinAmp)
                  << to_string(channel(i)->spinAmplitudes());

        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                std::cout << ",  ";
                std::static_pointer_cast<DecayingParticle>(d)->printDecayChainLevel(level + 1);
            }
    }

    if (level == 0)
        std::cout << "\n";
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
    auto bw = BlattWeisskopfs_.find(dt.freeAmplitude()->spinAmplitude()->L());
    if (bw == BlattWeisskopfs_.end())
        throw exceptions::Exception("No Blatt-Weisskopf factor found for L = "
                                    + std::to_string(dt.freeAmplitude()->spinAmplitude()->L()),
                                    "DecayingParticle::modifyDecayTree");

    // Add BlattWeisskopf object
    dt.addDataAccessor(bw->second.get());
}

//-------------------------
ComplexParameterVector DecayingParticle::freeAmplitudes() const
{
    ComplexParameterVector V;
    for (auto& c : Channels_) {
        auto vC = c->freeAmplitudes();
        V.insert(V.end(), vC.begin(), vC.end());
        // channels below
        for (auto& d : c->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                auto vD = std::static_pointer_cast<DecayingParticle>(d)->freeAmplitudes();
                V.insert(V.end(), vD.begin(), vD.end());
            }
    }

    // remove duplicates
    V.erase(ordered_unique(V.begin(), V.end()), V.end());

    return V;
}

//-------------------------
std::string DecayingParticle::printDecayTrees() const
{
    std::string s;
    for (const auto& m_dtv : DecayTrees_) {
        for (const auto& dt : m_dtv.second)
            s += "\ndepth = " + std::to_string(depth(*dt)) + "\n" + dt->asString() + "\n";
    }
    return s;
}

}
