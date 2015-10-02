#include "InitialStateParticle.h"

#include "Constants.h"
#include "DataSet.h"
#include "logging.h"

#include <TLorentzRotation.h>

#include <assert.h>

namespace yap {

//-------------------------
InitialStateParticle::InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    DecayingParticle(this, q, mass, name, radialSize),
    Prepared_(false),
    FourMomenta_(this),
    HelicityAngles_(this)
{
    addDataAccessor(this);

    // helicity angles do not store in Data_, so they don't need an index
    removeDataAccessor(&HelicityAngles_);
}

//-------------------------
InitialStateParticle::~InitialStateParticle()
{
    DataAccessors_.clear();
}

//-------------------------
double InitialStateParticle::logLikelihood(DataPartition& d)
{
    /// \todo implement

    DEBUG("InitialStateParticle::logLikelihood()");

    // update CalculationStatuses
    for (auto& pc : particleCombinations())
        updateCalculationStatus(d, pc);

    // calculate amplitues
    do {
        Amp a = Complex_0;
        for (auto& pc : particleCombinations()) {
            a += amplitude(d, pc);
        }

        DEBUG("InitialStateParticle amplitude = " << a);
    } while (d.increment());

    return 0;
}

//-------------------------
bool InitialStateParticle::consistent() const
{
    return DecayingParticle::consistent();
}

//-------------------------
bool InitialStateParticle::prepare()
{
    // check
    if (!consistent()) {
        LOG(ERROR) << "Cannot prepare InitialStateParticle, it is not consistent.";
        return false;
    }

    // particle combinations
    std::vector<std::shared_ptr<ParticleCombination> > PCs;
    for (auto& pc : particleCombinations()) {
        PCs.push_back(std::make_shared<ParticleCombination>(*pc));
    }
    ParticleCombination::makeParticleCombinationSetWithParents(PCs);
    clearSymmetrizationIndices();
    for (auto& pc : PCs)
        addSymmetrizationIndex(pc);

    // check
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (not pc->consistent()) {
            LOG(ERROR) << "Cannot prepare InitialStateParticle, particleCombinationSet is not consistent.";
            return false;
        }
        if (pc->indices().size() < particleCombinations()[0]->indices().size() and not pc->parent()) {
            LOG(ERROR) << "Cannot prepare InitialStateParticle, particleCombination is not consistent.";
            return false;
        }
    }

    //
    setSymmetrizationIndexParents();
    optimizeSpinAmplitudeSharing();

    // make sure that final state particles get the correct indices
    for (unsigned i = 0; i < particleCombinations()[0]->indices().size(); ++i) {
        std::shared_ptr<ParticleCombination> index = std::make_shared<ParticleCombination>(i);
        for (auto& pc : ParticleCombination::particleCombinationSet()) {
            if (ParticleCombination::equivByOrderlessContent(index, pc)) {
                FourMomenta_.addSymmetrizationIndex(pc);
                continue;
            }
        }
    }

    // add particle combinations to FourMomenta_ and HelicityAngles_
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        FourMomenta_.addSymmetrizationIndex(pc);

        if (pc->indices().size() > 1)
            HelicityAngles_.addSymmetrizationIndex(pc);
    }

    FourMomenta_.findInitialStateParticle();

    setDataAcessorIndices();

    // fill DecayChannels_
    DecayChannels_.clear();
    for (DataAccessor* d : DataAccessors_) {
        DecayChannel* ch = dynamic_cast<DecayChannel*>(d);
        if (ch)
            DecayChannels_.push_back(ch);
    }

    // check
    if (!consistent()) {
        LOG(ERROR) << "Something went wrong while preparing InitialStateParticle, it is not consistent anymore.";
        return false;
    }

    Prepared_ = true;

    return true;
}

//-------------------------
bool InitialStateParticle::setFreeAmplitudes(const std::vector<Amp>& amps)
{
    if (amps.size() != DecayChannels_.size()) {
        LOG(ERROR) << "InitialStateParticle::setFreeAmplitudes - amplitudes have wrong size "
                   << amps.size() << " != " << DecayChannels_.size();
        return false;
    }

    for (unsigned i = 0; i < amps.size(); ++i) {
        DecayChannels_[i]->setFreeAmplitude(amps[i]);
    }

    return true;
}

//-------------------------
std::vector<Amp> InitialStateParticle::freeAmplitudes() const
{
    std::vector<Amp> amps;
    amps.reserve(DecayChannels_.size());

    for (DecayChannel* ch : DecayChannels_)
        amps.push_back(ch->freeAmplitude());

    return amps;
}

//-------------------------
void InitialStateParticle::setSymmetrizationIndexParents()
{
    //std::cout << "InitialStateParticle::setSymmetrizationIndexParents()\n";

    unsigned size = particleCombinations()[0]->indices().size();
    assert(size > 1);

    clearSymmetrizationIndices();

    // get initial state PCs from set and add them
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (pc->indices().size() < size)
            continue;

        if (pc->daughters()[0]->parent() == nullptr)
            continue;

        //std::cout << "  add " << std::string(*pc) << "\n";
        addSymmetrizationIndex(pc);
    }

    // next level
    for (auto& ch : channels())
        ch->setSymmetrizationIndexParents();

}

//-------------------------
bool InitialStateParticle::addDataPoint(const std::vector<TLorentzVector>& fourMomenta)
{
    if (!Prepared_) {
        LOG(ERROR) << "Cannot add DataPoint to InitialStateParticle. Call InitialStateParticle::prepare() first!";
        return false;
    }

    if (DataSet_.empty()) {
        return addDataPoint(DataPoint(fourMomenta));
    }

    DataSet_.push_back(DataPoint(DataSet_[0]));

    DataPoint& d = DataSet_.back();

    if (! d.setFourMomenta(fourMomenta))
        return false;

    FourMomenta_.calculate(d);
    HelicityAngles_.calculate(d);

    if (!DataSet_.consistent(d))
        return false;

    return true;
}

//-------------------------
bool InitialStateParticle::addDataPoint(DataPoint&& d)
{
    if (!Prepared_) {
        LOG(ERROR) << "Cannot add DataPoint to InitialStateParticle. Call InitialStateParticle::prepare() first!";
        return false;
    }

    d.allocateStorage(FourMomenta_, HelicityAngles_, DataAccessors_);

    FourMomenta_.calculate(d);
    HelicityAngles_.calculate(d);

    /*for (auto& pc : FourMomenta_.particleCombinations()) {
      std::cout << FourMomenta_.symmetrizationIndex(pc) << "  \t" << std::string(*pc) << "  \t";
      FourMomenta_.p(d, FourMomenta_.symmetrizationIndex(pc)).Print();
    }

    for (auto& pc : HelicityAngles_.particleCombinations()) {
      std::cout << HelicityAngles_.symmetrizationIndex(pc) << "  \t" << std::string(*pc) << "  \t";
      for (double a : HelicityAngles_.helicityAngles(d, HelicityAngles_.symmetrizationIndex(pc)))
          std::cout << a << "  \t";
      std::cout << "\n";
    }*/

    d.printDataSize();

    if (!DataSet_.consistent(d))
        return false;

    DataSet_.push_back(d);
    return true;
}

//-------------------------
bool InitialStateParticle::addDataPoint(const DataPoint& d)
{
    return addDataPoint(DataPoint(d));
}

//-------------------------
void InitialStateParticle::printDataAccessors()
{
    std::cout << "DataAccessors of " << name() << "\n"
              << "index \tnSymIndices \taddress  \tname  \t\tparticleCombinations\n";
    for (DataAccessor* d : DataAccessors_) {
        std::cout << d->index() << "  \t" << d->maxSymmetrizationIndex() + 1
                  << "  \t" << d
                  << "  \t(" << typeid(*d).name() << ")  \t";
        if (dynamic_cast<Particle*>(d))
            std::cout << dynamic_cast<Particle*>(d)->name();
        else if (dynamic_cast<DecayChannel*>(d))
            std::cout << std::string(*dynamic_cast<DecayChannel*>(d));

        std::cout << " \t";

        for (auto& pc : d->particleCombinations())
            std::cout << std::string(*pc) << ";  ";
        std::cout << "\n";
    }
    std::cout << std::endl;
}

//-------------------------
void InitialStateParticle::removeDataAccessor(DataAccessor* d)
{
    if (! DataAccessors_.empty())
        DataAccessors_.erase(d);

    //DEBUG("size of InitialStateParticle's DataAccessors_ = " << DataAccessors_.size());
}

//-------------------------
void InitialStateParticle::setDataAcessorIndices()
{
    unsigned i(0);
    for (DataAccessor* d : DataAccessors_)
        d->setIndex(i++);
}


}
