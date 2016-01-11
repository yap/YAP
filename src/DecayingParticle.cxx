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
    DataAccessor(),
    RadialSize_(std::make_shared<RealParameter>(radialSize)),
    Amplitude_(std::make_shared<ComplexCachedDataValue>(this))
{
}

//-------------------------
std::complex<double> DecayingParticle::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex) const
{
    unsigned symIndex = symmetrizationIndex(pc);

    if (Amplitude_->calculationStatus(pc, symIndex, dataPartitionIndex) == kUncalculated) {

        std::complex<double> a = Complex_0;

        /// \todo Is this the best way to do it? (loop over pc's then channels. or channels then pc's?)

        // sum up DecayChannel::amplitude over each channel
        for (auto& c : channels()) {
            if (c->hasSymmetrizationIndex(pc))
                a += c->amplitude(d, pc, dataPartitionIndex);
        }

        Amplitude_->setValue(a, d, symIndex, dataPartitionIndex);

        DEBUG("DecayingParticle::amplitude - calculate amplitude for " << name() << " " << *pc << " = " << a);
        return a;
    }

    DEBUG("DecayingParticle::amplitude - using cached amplitude for " << name() << " " << *pc << " = " << Amplitude_->value(d, symIndex));
    return Amplitude_->value(d, symIndex);
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
            BlattWeisskopfs_.insert(std::make_pair(L, std::make_shared<BlattWeisskopf>(L, this)));

        Channels_.back()->amplitudes(sa).Fixed->addDependencies(BlattWeisskopfs_[sa->L()]->ParametersItDependsOn());
        Channels_.back()->amplitudes(sa).Fixed->addDependencies(BlattWeisskopfs_[sa->L()]->CachedDataValuesItDependsOn());
    }
            
    
    // add particle combinations
    for (auto pc : Channels_.back()->particleCombinations()) {
        addSymmetrizationIndex(pc);
        // and add pc's daughters to Channels' daughters
        // for (size_t i = 0; i < pc->daughters().size(); ++i)
        //     Channels_.back()->Daughters_[i].addSymmetrizationIndex(pc->daughters()[i]);
    }

    // add dependencies
    Amplitude_->addDependencies(Channels_.back()->ParametersItDependsOn());
    Amplitude_->addDependencies(Channels_.back()->CachedDataValuesItDependsOn());

    FLOG(INFO) << *Channels_.back() << " with N(PC) = " << Channels_.back()->particleCombinations().size();
}

//-------------------------
void DecayingParticle::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> pc)
{
    DataAccessor::addSymmetrizationIndex(pc);
    // add also to all BlattWeiskopf barrier factors
    for (auto& kv : BlattWeisskopfs_)
        kv.second->addSymmetrizationIndex(pc);
}

//-------------------------
void DecayingParticle::addChannel(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, unsigned l)
{
    // check isp's are set
    if (!A->initialStateParticle() or !B->initialStateParticle())
        throw exceptions::Exception("A daughter has an unset InitialStateParticle", "DecayingParticle::addChannel");
    // check isp's are the same
    if (A->initialStateParticle() != B->initialStateParticle())
        throw exceptions::Exception("InitialStateParticle mismatch among daughters", "DecayingParticle::addChannel");
    // and check isp of this is either unset or is the same as the daughters
    if (initialStateParticle() and A->initialStateParticle() != initialStateParticle())
        throw exceptions::Exception("InitialStateParticle mismatch of daughters", "DecayingParticle::addChannel");

    auto sa = A->initialStateParticle()->spinAmplitudeCache().spinAmplitude(quantumNumbers(), A->quantumNumbers(), B->quantumNumbers(), l);
    addChannel(std::make_unique<DecayChannel>(ParticleVector{A, B}, sa));
}

//-------------------------
void DecayingParticle::addChannels(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, unsigned max_l)
{
    for (unsigned l = 0; l <= max_l; ++l) {
        try { addChannel(A, B, l);}
        catch (const exceptions::AngularMomentumNotConserved&) {/* ignore */}
        catch (const exceptions::InconsistentSpinProjection&) {/* ignore */}
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
                paddingSpinAmp = std::max(paddingSpinAmp, std::string(*c->spinAmplitude()).length());
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
                  << std::string(*(channel(i)->spinAmplitude()));

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
void DecayingParticle::printSpinAmplitudes(int level)
{
    static std::set<std::shared_ptr<SpinAmplitude> > ampSet;

    // get maximum length of particle names
    static size_t padding = 0;
    static size_t paddingSpinAmp = 6;
    if (padding == 0 || level == -1) {
        ampSet.clear();
        padding = std::max(padding, name().length());
        for (auto& c : Channels_) {
            for (std::shared_ptr<Particle> d : c->daughters()) {
                padding = std::max(padding, d->name().length());
                if (std::dynamic_pointer_cast<DecayingParticle>(d))
                    std::static_pointer_cast<DecayingParticle>(d)->printSpinAmplitudes(-1);
                paddingSpinAmp = std::max(paddingSpinAmp, std::string(*c->spinAmplitude()).length());
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

        std::shared_ptr<SpinAmplitude> amp = channel(i)->spinAmplitude();
        if (ampSet.find(amp) == ampSet.end()) {
            ampSet.insert(amp);
            std::cout << std::left << std::setw(paddingSpinAmp) << std::string(*amp);
        } else {
            std::cout << std::left << std::setw(paddingSpinAmp - 1) << "shared";
        }


        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                std::cout << ",  ";
                std::static_pointer_cast<DecayingParticle>(d)->printSpinAmplitudes(level + 1);
            }
    }

    if (level == 0)
        std::cout << "\n";
}

//-------------------------
void DecayingParticle::setSymmetrizationIndexParents()
{
    // clean up PCs without parents
    auto PCsParents = particleCombinations();
    auto it = PCsParents.begin();
    while (it != PCsParents.end()) {
        if (!(*it)->parent()) {
            it = PCsParents.erase(it);
        } else
            ++it;
    }
    clearSymmetrizationIndices();

    for (auto& pc : PCsParents)
        addSymmetrizationIndex(pc);

    // next level
    for (auto& ch : channels())
        ch->setSymmetrizationIndexParents();

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
        V.emplace_back(c->freeAmplitude());
        // add channels below
        for (auto& d : c->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                auto v = std::dynamic_pointer_cast<DecayingParticle>(d)->freeAmplitudes();
                V.insert(V.end(), v.begin(), v.end());
            }
    }

    // remove duplicates
    V.erase(ordered_unique(V.begin(), V.end()), V.end());

    return V;
}

}
