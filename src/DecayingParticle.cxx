#include "DecayingParticle.h"

#include "container_utils.h"
#include "HelicitySpinAmplitude.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "ParticleCombinationCache.h"

#include <iomanip>
#include <memory>
#include <stdexcept>

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
std::complex<double> DecayingParticle::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, int two_m, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    // get cached amplitude object for spin projection two_m
    auto A = Amplitudes_.at(two_m);

    if (A->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        std::complex<double> a = Complex_0;

        /// \todo Is this the best way to do it? (loop over pc's then channels. or channels then pc's?)

        // sum up DecayChannel::amplitude over each channel
        for (auto& c : channels()) {
            if (c->hasParticleCombination(pc))
                a += c->amplitude(d, pc, two_m, dataPartitionIndex);
        }

        A->setValue(a, d, symIndex, dataPartitionIndex);

        DEBUG("DecayingParticle::amplitude - calculate amplitude for " << name() << " " << *pc << " = " << a);
        return a;
    }

    DEBUG("DecayingParticle::amplitude - using cached amplitude for " << name() << " " << *pc << " = " << A->value(d, symIndex));
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
void DecayingParticle::addChannel(std::unique_ptr<DecayChannel> c)
{
    if (!c)
        throw exceptions::Exception("DecayChannel empty", "DecayingParticle::addChannel");

    if (c->particleCombinations().empty())
        throw exceptions::Exception(std::string("DecayChannel has no ParticleCombinations - ") + to_string(*c),
                                    "DecayingParticle::addChannel");

    // check ISP
    if (!Channels_.empty() and c->initialStateParticle() != initialStateParticle())
        throw exceptions::Exception("InitialStateParticle mismatch", "DecayingParticle::addChannel");

    Channels_.emplace_back(std::move(c));
    Channels_.back()->setDecayingParticle(this);

    // insert necessary Blatt-Weisskopf barrier factor
    // and set dependencies for DecayChannel dependent on it
    for (auto& sa : Channels_.back()->spinAmplitudes()) {

        // if BW is not already stored for L, add it
        if (BlattWeisskopfs_.find(sa->L()) == BlattWeisskopfs_.end())
            BlattWeisskopfs_.insert(std::make_pair(sa->L(), std::make_shared<BlattWeisskopf>(sa->L(), this)));

        // add BW to Fixed amplitudes for all spin projections
        auto& apM = Channels_.back()->amplitudes(sa);
        for (auto& ap : apM) {
            ap.second.Fixed->addDependencies(BlattWeisskopfs_[sa->L()]->ParametersItDependsOn());
            ap.second.Fixed->addDependencies(BlattWeisskopfs_[sa->L()]->CachedDataValuesItDependsOn());
        }
    }

    // add particle combinations
    for (auto pc : Channels_.back()->particleCombinations()) {
        addParticleCombination(pc);
        // and add pc's daughters to Channels' daughters
        // for (size_t i = 0; i < pc->daughters().size(); ++i)
        //     Channels_.back()->Daughters_[i].addSymmetrizationIndex(pc->daughters()[i]);
    }

    // Add DecayChannel's TotalAmplitude's as dependencies for this object's Amplitudes
    // by spin projection (key in TotalAmplitudes_)
    for (auto& kv : Channels_.back()->TotalAmplitudes_) {
        // if spin projection not yet in Amplitudes_, add it
        if (Amplitudes_.find(kv.first) == Amplitudes_.end())
            Amplitudes_[kv.first] = ComplexCachedDataValue::create(this);
        Amplitudes_[kv.first]->addDependency(kv.second);
    }

    FLOG(INFO) << *Channels_.back() << " with N(PC) = " << Channels_.back()->particleCombinations().size();
}

//-------------------------
void DecayingParticle::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    DataAccessor::addParticleCombination(pc);

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
}

//-------------------------
std::vector< std::shared_ptr<FinalStateParticle> > DecayingParticle::finalStateParticles(unsigned i) const
{
    if (!Channels_.at(i))
        return std::vector<std::shared_ptr<FinalStateParticle>>();
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
CachedDataValueSet DecayingParticle::CachedDataValuesItDependsOn()
{
    CachedDataValueSet S;
    for (auto& kv : Amplitudes_)
        S.insert(kv.second);
    return S;
}

//-------------------------
DataAccessorSet DecayingParticle::dataAccessors()
{
    DataAccessorSet V;
    for (auto& c : Channels_) {
        // add channel
        V.emplace(c);
        // and channel's data accessors
        auto v = c->dataAccessors();
        V.insert(v.begin(), v.end());
    }
    return V;
}

//-------------------------
ComplexParameterVector DecayingParticle::freeAmplitudes() const
{
    ComplexParameterVector V;
    for (auto& c : Channels_) {
        // add channel
        auto vC = c->freeAmplitudes();
        V.insert(V.end(), vC.begin(), vC.end());
        // add channels below
        for (auto& d : c->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                auto vD = std::dynamic_pointer_cast<DecayingParticle>(d)->freeAmplitudes();
                V.insert(V.end(), vD.begin(), vD.end());
            }
    }

    // remove duplicates
    V.erase(ordered_unique(V.begin(), V.end()), V.end());

    return V;
}

}
