#include "InitialStateParticle.h"

#include "Constants.h"
#include "DataPartition.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "logging.h"

#include <assert.h>
#include <future>

namespace yap {

//-------------------------
InitialStateParticle::InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    DecayingParticle(q, mass, name, radialSize),
    Prepared_(false),
    CoordinateSystem_(ThreeAxes),
    FourMomenta_(),
    MeasuredBreakupMomenta_(),
    HelicityAngles_()
{
    FourMomenta_.setInitialStateParticle(this);
    MeasuredBreakupMomenta_.setInitialStateParticle(this);
    HelicityAngles_.setInitialStateParticle(this);

    setInitialStateParticle(this);
}

//-------------------------
InitialStateParticle::~InitialStateParticle()
{
    DataAccessors_.clear();
}

//-------------------------
std::complex<double> InitialStateParticle::amplitude(DataPoint& d, unsigned dataPartitionIndex) const
{
    std::complex<double> a = Complex_0;
    for (auto& pc : particleCombinations())
        a += amplitude(d, pc, dataPartitionIndex);

    DEBUG ("InitialStateParticle::amplitude = " << a);

    return a;
}

//-------------------------
double InitialStateParticle::partialSumOfLogsOfSquaredAmplitudes(DataPartitionBase* D)
{
    if (!hasDataPartition(D))
        return 0;

    double L = 0;

    // loop over data points in partition
    for (DataIterator d = D->begin(); d != D->end(); ++d) {
        DEBUG("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ");

        L += logOfSquaredAmplitude(*d, D->index());
    }

    setCachedDataValueFlagsToUnchanged(D->index());

    return L;
}

//-------------------------
double InitialStateParticle::sumOfLogsOfSquaredAmplitudes()
{
    // update global caclulationStatus's before looping over partitions
    updateGlobalCalculationStatuses();

    // create thread for calculation on each partition
    std::vector<std::future<double> > L;
    for (auto& D : DataPartitions_)
        L.push_back(std::async(std::launch::async, &InitialStateParticle::partialSumOfLogsOfSquaredAmplitudes, this, D.get()));

    // wait for each partition to finish calculating
    double sum = 0;
    for (auto& l : L)
        sum += l.get();

    setParameterFlagsToUnchanged();

    return sum;
}

//-------------------------
bool InitialStateParticle::consistent() const
{
    bool result(true);

    result &= DecayingParticle::consistent();
    result &= FourMomenta_.consistent();
    result &= MeasuredBreakupMomenta_.consistent();
    result &= HelicityAngles_.consistent();

    /// \todo: is this necessary to check?
    /// @{
    std::vector<std::shared_ptr<FinalStateParticle> > A = FinalStateParticles_;
    sort(A.begin(), A.end());
    std::vector<std::shared_ptr<FinalStateParticle> > B = DecayingParticle::finalStateParticles();
    sort(B.begin(), B.end());

    if (A != B) {
        FLOG(ERROR) << "FinalStateParticles_ are not set correctly.";
        result = false;
    }
    /// @}

    return result;
}

//-------------------------
bool InitialStateParticle::prepare()
{
    // check
    if (!DecayingParticle::consistent()) {
        FLOG(ERROR) << "Cannot prepare InitialStateParticle, it is not consistent.";
        return false;
    }

    // check coordinate system
    CoordinateSystem_ = unit(CoordinateSystem_);
    if (!isRightHanded(CoordinateSystem_)) {
        FLOG(ERROR) << "Coordinate system is not right-handed.";
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

    // add non-final-state particle combinations to FourMomenta_, MeasuredBreakupMomenta_ and HelicityAngles_
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (pc->isFinalStateParticle())
            continue;
        FourMomenta_.addSymmetrizationIndex(pc);
        HelicityAngles_.addSymmetrizationIndex(pc);
        MeasuredBreakupMomenta_.addSymmetrizationIndex(pc);
    }

    // prepare FourMomenta. Needs FinalStateParticles_
    FourMomenta_.prepare();

    // set consecutive indices for DataAccessors_
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
bool InitialStateParticle::setFinalStateParticles(std::initializer_list<std::shared_ptr<FinalStateParticle> > FSP)
{
    // check that FinalStateParticles_ is empty
    if (!FinalStateParticles_.empty()) {
        FLOG(ERROR) << "final-state particles have already been set.";
        return false;
    }

    // check that none of the FSP's has yet been used
    bool all_unused = true;
    for (auto& fsp : FSP)
        if (!fsp->particleCombinations().empty()) {
            FLOG(ERROR) << "final-state particle already has indices assigned: " << (std::string)*fsp;
            all_unused = false;
        }
    if (!all_unused)
        return false;

    FinalStateParticles_.reserve(FSP.size());

    // set indices by order in vector
    for (auto& fsp : FSP) {
        fsp->addSymmetrizationIndex(ParticleCombination::uniqueSharedPtr(FinalStateParticles_.size()));
        FinalStateParticles_.push_back(fsp);
    }

    return true;
}

//-------------------------
std::array<double, 2> InitialStateParticle::getMassRange(const std::shared_ptr<const ParticleCombination>& pc) const
{
    std::array<double, 2> m = {0, mass()->value()};

    for (size_t i = 0; i < FinalStateParticles_.size(); ++i) {
        if (std::find(pc->indices().begin(), pc->indices().end(), i) != pc->indices().end())
            // add mass to low end
            m[0] += FinalStateParticles_[i]->mass()->value();
        else
            // subtract mass from high end
            m[1] -= FinalStateParticles_[i]->mass()->value();
    }
    return m;
}

//-------------------------
std::vector<std::shared_ptr<ComplexParameter> > InitialStateParticle::freeAmplitudes() const
{
    std::vector<std::shared_ptr<ComplexParameter> > amps;
    amps.reserve(DecayChannels_.size());

    for (DecayChannel* ch : DecayChannels_) {
        // if there is only one channel, its amplitude is 1
        // no need to change/fit it
        if (ch->parent()->nChannels() == 1)
            continue;

        amps.push_back(ch->freeAmplitude());
    }

    return amps;
}

//-------------------------
std::vector<DataPartitionBase*> InitialStateParticle::dataPartitions()
{
    std::vector<DataPartitionBase*> partitions;
    partitions.reserve(DataPartitions_.size());

    for (auto& d : DataPartitions_)
        partitions.push_back(d.get());

    return partitions;
}

//-------------------------
void InitialStateParticle::setDataPartitions(std::vector<std::unique_ptr<DataPartitionBase> > partitions)
{
    DataPartitions_ = std::move(partitions);

    for (unsigned i = 0; i < DataPartitions_.size(); ++i) {
        DataPartitions_[i]->setIndex(i);
    }

    setNumberOfDataPartitions(DataPartitions_.size());
}

//-------------------------
bool InitialStateParticle::addDataPoint(const std::vector<FourVector<double> >& fourMomenta)
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

    if (! d.setFinalStateFourMomenta(fourMomenta))
        return false;

    /// calculate FourMomenta_, MeasuredBreakupMomenta_ and HelicityAngles_
    calculate(d);

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

    d.allocateStorage(FourMomenta_, DataAccessors_);

    /// calculate FourMomenta_, MeasuredBreakupMomenta_ and HelicityAngles_
    calculate(d);

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
bool InitialStateParticle::initializeForMonteCarloGeneration(unsigned n)
{
    bool result = true;

    // initialize with 0
    std::vector<FourVector<double> > momenta(FinalStateParticles_.size(), {0.,0.,0.,0.});

    // add n (empty) data points
    for (unsigned i = 0; i < n; ++i)
        result &= addDataPoint(momenta);

    // set data partitions (1 for each data point)
    setDataPartitions(createDataPartitionsBlockBySize(DataSet_, 1));

    // do one initial calculation
    /// \todo Only calculate data-independent values
    result &= std::isfinite(sumOfLogsOfSquaredAmplitudes());

    return result;
}

//-------------------------
void InitialStateParticle::printDataAccessors(bool printParticleCombinations)
{
    // header
    std::cout << "DataAccessors of " << name() << "\n"
              << "index \tnSymIndices \taddress  \tname";
    if (printParticleCombinations)
        std::cout << "\t\tparticleCombinations";
    std::cout << "\n";

    for (DataAccessor* d : DataAccessors_) {
        std::cout << d->index() << "  \t" << d->maxSymmetrizationIndex() + 1
                  << "  \t\t" << d
                  << "  \t(" << typeid(*d).name() << ")  \t";
        if (dynamic_cast<Particle*>(d))
            std::cout << dynamic_cast<Particle*>(d)->name();
        else if (dynamic_cast<DecayChannel*>(d))
            std::cout << std::string(*dynamic_cast<DecayChannel*>(d));

        if (printParticleCombinations) {
            std::cout << " \t";

            for (auto& pc : d->particleCombinations())
                std::cout << std::string(*pc) << ":" << d->symmetrizationIndex(pc) << ";  ";
        }

        std::cout << "\n";
    }
    std::cout << std::endl;
}

//-------------------------
void InitialStateParticle::setSymmetrizationIndexParents()
{
    unsigned size = particleCombinations()[0]->indices().size();
    assert(size > 1);

    clearSymmetrizationIndices();

    // get initial state PCs from set and add them
    for (auto& pc : ParticleCombination::particleCombinationSet()) {
        if (pc->indices().size() < size)
            continue;

        if (pc->daughters()[0]->parent() == nullptr)
            continue;

        addSymmetrizationIndex(pc);
    }

    // next level
    for (auto& ch : channels())
        ch->setSymmetrizationIndexParents();

}

//-------------------------
bool InitialStateParticle::hasDataPartition(DataPartitionBase* d)
{
    for (auto& dp : DataPartitions_) {
        if (d == dp.get())
            return true;
    }

    LOG(ERROR) << "InitialStateParticle::hasDataPartition - trying to calculate for a DataPartition, which is not stored in this InitialStateParticle.";
    return false;
}

//-------------------------
void InitialStateParticle::setNumberOfDataPartitions(unsigned n)
{
    for (DataAccessor* d : DataAccessors_) {
        if (d == &FourMomenta_ or d == &MeasuredBreakupMomenta_  or d == &HelicityAngles_)
            continue;

        for (CachedDataValue* c : d->CachedDataValues_) {
            c->setNumberOfDataPartitions(n);
        }
    }
}

//-------------------------
void InitialStateParticle::updateGlobalCalculationStatuses()
{
    for (DataAccessor* d : DataAccessors_) {
        /// \todo Think hard about a less stupid solution
        if (d == this)
            continue;
        d->updateGlobalCalculationStatuses();
    }
}

//-------------------------
void InitialStateParticle::calculate(DataPoint& d)
{
    FourMomenta_.calculate(d);
    MeasuredBreakupMomenta_.calculate(d);
    HelicityAngles_.calculate(d);
}

//-------------------------
void InitialStateParticle::resetCalculationStatuses(unsigned dataPartitionIndex)
{
    for (DataAccessor* d : DataAccessors_) {
        if (d == &FourMomenta_ or d == &MeasuredBreakupMomenta_  or d == &HelicityAngles_)
            continue;

        //DEBUG("resetCalculationStatus for " << typeid(*d).name() << "  " << d);

        for (CachedDataValue* c : d->CachedDataValues_) {
            c->resetCalculationStatus(dataPartitionIndex);
        }
    }
}

//-------------------------
void InitialStateParticle::setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex)
{
    for (DataAccessor* d : DataAccessors_) {
        if (d == &FourMomenta_ or d == &MeasuredBreakupMomenta_  or d == &HelicityAngles_)
            continue;

        for (CachedDataValue* c : d->CachedDataValues_) {
            c->setVariableStatus(kUnchanged, dataPartitionIndex);
        }
    }
}

//-------------------------
void InitialStateParticle::setParameterFlagsToUnchanged()
{
    for (DataAccessor* d : DataAccessors_) {
        if (d == &FourMomenta_ or d == &MeasuredBreakupMomenta_  or d == &HelicityAngles_)
            continue;

        for (CachedDataValue* c : d->CachedDataValues_) {
            for (auto& p : c->ParametersItDependsOn_) {
                if (p->variableStatus() == kChanged) {
                    p->setVariableStatus(kUnchanged);
                }
            }
        }
    }
}

//-------------------------
void InitialStateParticle::addDataAccessor(DataAccessor* d)
{
    DataAccessors_.insert(d);
    if (prepared()) {
        LOG(ERROR) << "InitialStateParticle has already been prepared. "
                   << "Do NOT modify/add DecayChannels etc. after calling InitialStateParticle::prepare(), "
                   << "otherwise it will become inconsistent!";
    }
}


//-------------------------
void InitialStateParticle::removeDataAccessor(DataAccessor* d)
{
    if (! DataAccessors_.empty()) {

        DataAccessors_.erase(d);

        if (prepared()) {
            LOG(ERROR) << "InitialStateParticle has already been prepared. "
                       << "Do NOT modify/add DecayChannels etc. after calling InitialStateParticle::prepare(), "
                       << "otherwise it will become inconsistent!";
        }
    }
}

//-------------------------
void InitialStateParticle::setDataAcessorIndices()
{
    unsigned i(0);
    for (DataAccessor* d : DataAccessors_)
        d->setIndex(i++);
}


}
