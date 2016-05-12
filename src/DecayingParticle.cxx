#include "DecayingParticle.h"

#include "AmplitudePair.h"
#include "BlattWeisskopf.h"
#include "CalculationStatus.h"
#include "Constants.h"
#include "container_utils.h"
#include "DecayChannel.h"
#include "DecayTree.h"
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
    AmplitudeComponent(),
    Particle(q, mass, name),
    DataAccessor(&ParticleCombination::equivUpAndDown),
    RadialSize_(std::make_shared<RealParameter>(radialSize))
{
}

//-------------------------
std::complex<double> DecayingParticle::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, int two_m, StatusManager& sm) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    // get cached amplitude object for spin projection two_m
    auto& A = Amplitudes_.at(two_m);

    if (sm.status(*A, symIndex) == CalculationStatus::uncalculated) {

        std::complex<double> a = Complex_0;

        /// \todo Is this the best way to do it? (loop over pc's then channels. or channels then pc's?)

        // sum up DecayChannel::amplitude over each channel
        for (const auto& c : channels())
            if (c->hasParticleCombination(pc))
                a += c->amplitude(d, pc, two_m, sm);

        A->setValue(a, d, symIndex, sm);

        DEBUG("DecayingParticle::amplitude - calculated amplitude for " << name() << " " << *pc << " = " << a);
        return a;
    }

    DEBUG("DecayingParticle::amplitude - used cached amplitude for " << name() << " " << *pc << " = " << A->value(d, symIndex));
    return A->value(d, symIndex);
}

//-------------------------
bool DecayingParticle::consistent() const
{
    bool C = DataAccessor::consistent();
    C &= Particle::consistent();

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
    addToModel();

    // if this is to be the initial state particle
    if (!model()->initialStateParticle() and finalStateParticles().size() == model()->finalStateParticles().size())
        const_cast<Model*>(static_cast<const DecayingParticle*>(this)->model())->setInitialStateParticle(std::static_pointer_cast<DecayingParticle>(shared_from_this()));

    // add particle combinations
    for (auto pc : Channels_.back()->particleCombinations()) {
        addParticleCombination(pc);
    }

    // Add DecayChannel's TotalAmplitude's as dependencies for this object's Amplitudes
    // by spin projection (key in TotalAmplitudes_)
    for (auto& kv : Channels_.back()->TotalAmplitudes_) {
        // if spin projection not yet in Amplitudes_, add it
        if (Amplitudes_.find(kv.first) == Amplitudes_.end())
            Amplitudes_[kv.first] = ComplexCachedDataValue::create(this);
        Amplitudes_[kv.first]->addDependency(kv.second);
    }

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
unsigned DecayingParticle::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    unsigned index = DataAccessor::addParticleCombination(pc);

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

    return index;
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
CachedDataValueSet DecayingParticle::cachedDataValuesItDependsOn()
{
    CachedDataValueSet S;
    for (auto& kv : Amplitudes_)
        S.insert(kv.second);
    return S;
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
DecayTreeVector DecayingParticle::decayTrees(int two_m) const
{
    DecayTreeVector dtv;

    // loop over decay channels
    for (const auto& c : Channels_) {
        // loop over entries in map (spin amplitude -> decay tree vector) for spin projection
        for (const auto& sa_dtv : c->DecayTrees_[two_m]) {
            // insert DecayTree's
            dtv.insert(dtv.end(), sa_dtv.second.begin(), sa_dtv.second.end());
        }
    }

    return dtv;
}

//-------------------------
void DecayingParticle::modifyDecayTree(std::shared_ptr<DecayTree> dt) const
{
    if (!dt)
        throw exceptions::Exception("DecayTree is nullptr", "DecayingParticle::modifyDecayTree");

    // find SpinAmplitude in DecayTree
    const SpinAmplitude* sa = nullptr;
    for (const auto& sda : dt->StaticDataAccessors_)
        if (dynamic_cast<const SpinAmplitude*>(sda)) {
            sa = static_cast<const SpinAmplitude*>(sda);
            break;
        }
    if (!sa)
        throw exceptions::Exception("DecayTree has no SpinAmplitude", "DecayingParticle::modifyDecayTree");

    // find BlattWeisskopf object
    auto bw = BlattWeisskopfs_.find(sa->L());
    if (bw == BlattWeisskopfs_.end())
        throw exceptions::Exception("No Blatt-Weisskopf factor found for L = " + std::to_string(sa->L()),
                                    "DecayingParticle::modifyDecayTree");

    // Add BlattWeisskopf object
    dt->addDataAccessor(bw->second.get());
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

}
