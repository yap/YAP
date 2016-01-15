#include "InitialStateParticle.h"

#include "Constants.h"
#include "DataPartition.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "ParticleCombinationCache.h"

#include <assert.h>
#include <future>
#include <stdexcept>

namespace yap {

//-------------------------
InitialStateParticle::InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    std::enable_shared_from_this<InitialStateParticle>(),
    DecayingParticle(q, mass, name, radialSize),
    Prepared_(false),
    CoordinateSystem_(ThreeAxes),
    SpinAmplitudeCache_(this),
    FourMomenta_(std::make_shared<FourMomenta>(this)),
    MeasuredBreakupMomenta_(std::make_shared<MeasuredBreakupMomenta>(this)),
    HelicityAngles_(std::make_shared<HelicityAngles>(this))
{
}

//-------------------------
std::complex<double> InitialStateParticle::amplitude(DataPoint& d, int two_m, unsigned dataPartitionIndex) const
{
    std::complex<double> a = Complex_0;

    // sum up DecayingParticle::amplitude over each particle combination
    for (auto& pc : particleCombinations())
        a += amplitude(d, pc, two_m, dataPartitionIndex);

    return a;
}

//-------------------------
std::complex<double> InitialStateParticle::amplitude(DataPoint& d, unsigned dataPartitionIndex) const
{
    std::complex<double> a = Complex_0;

    for (int two_m = -quantumNumbers().twoJ(); two_m <= (int)quantumNumbers().twoJ(); two_m += 2)
        a += amplitude(d, two_m, dataPartitionIndex);

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
    bool C = DecayingParticle::consistent();
    C &= FourMomenta_->consistent();
    C &= MeasuredBreakupMomenta_->consistent();
    C &= HelicityAngles_->consistent();
    C &= ParticleCombinationCache_.consistent();
    C &= SpinAmplitudeCache_.consistent();

    /// \todo: is this necessary to check?
    /// @{
    std::vector<std::shared_ptr<FinalStateParticle> > A = FinalStateParticles_;
    sort(A.begin(), A.end());
    std::vector<std::shared_ptr<FinalStateParticle> > B = DecayingParticle::finalStateParticles();
    sort(B.begin(), B.end());

    if (A != B) {
        FLOG(ERROR) << "FinalStateParticles_ are not set correctly.";
        C &= false;
    }
    /// @}

    return C;
}

//-------------------------
void InitialStateParticle::prepare()
{
    // check
    if (!DecayingParticle::consistent()) {
        FLOG(ERROR) << "Cannot prepare InitialStateParticle, it is not consistent as DecayingParticle.";
        throw exceptions::Exception("InitialStateParticle inconsistent", "InitialStateParticle::prepare");
    }

    // check coordinate system
    CoordinateSystem_ = unit(CoordinateSystem_);
    if (!isRightHanded(CoordinateSystem_)) {
        FLOG(ERROR) << "Coordinate system is not right-handed.";
        throw exceptions::Exception("Coordinate system not right-handed", "InitialStateParticle::prepare");
    }

    /*
        // particle combinations
        std::vector<std::shared_ptr<ParticleCombination> > PCs;
        for (auto& pc : particleCombinations()) {
            PCs.push_back(std::make_shared<const ParticleCombination>(*pc));
        }
        particleCombinationCache = ParticleCombinationCache(PCs);
        clearSymmetrizationIndices();
        for (auto& pc : PCs)
            addSymmetrizationIndex(pc);
    */
    // check
    for (auto& wpc : ParticleCombinationCache_) {
        if (!wpc.lock())
            continue;
        auto pc = wpc.lock();
        if (!pc->consistent()
                or (pc->indices().size() < particleCombinations()[0]->indices().size() and !pc->parent())) {
            FLOG(ERROR) << "Cannot prepare InitialStateParticle, particleCombinationCache is not consistent.";
            throw exceptions::Exception("ParticleCombination inconsistent", "InitialStateParticle::prepare");
        }
    }

    //
    setSymmetrizationIndexParents();

    // add non-final-state particle combinations to FourMomenta_, MeasuredBreakupMomenta_ and HelicityAngles_
    for (auto& wpc : ParticleCombinationCache_) {
        if (wpc.expired())
            continue;
        auto pc = wpc.lock();
        if (pc->isFinalStateParticle())
            continue;
        FourMomenta_->addSymmetrizationIndex(pc);
        HelicityAngles_->addSymmetrizationIndex(pc);
        MeasuredBreakupMomenta_->addSymmetrizationIndex(pc);
    }

    // prepare FourMomenta. Needs FinalStateParticles_
    FourMomenta_->prepare();

    // set DataAccessors_
    DataAccessors_ = dataAccessors();
    // add this (commented out because ISP has no need for data access at moment)
    // DataAccessors_.push_back(shared_from_this());
    // set unique indices to all DataAccessors
    unsigned i = 0;
    for (auto da : DataAccessors_)
        da->setIndex(i++);

    // check
    if (!consistent()) {
        FLOG(ERROR) << "Something went wrong while preparing InitialStateParticle, it is not consistent anymore.";
        throw exceptions::Exception("InitialStateParticle inconsistent", "InitialStateParticle::prepare");
    }

    Prepared_ = true;
}

//-------------------------
void InitialStateParticle::setFinalStateParticles(std::initializer_list<std::shared_ptr<FinalStateParticle> > FSP)
{
    // check that FinalStateParticles_ is empty
    if (!FinalStateParticles_.empty()) {
        FLOG(ERROR) << "final-state particles have already been set.";
        throw exceptions::Exception("Final-state particles already set", "InitialStateParticle::setFinalStateParticles");
    }

    // check that none of the FSP's has yet been used
    // and that FinalStateParticles don't have ISP set to another ISP
    // (in principle ISP's should be set to nullptr)
    for (auto& fsp : FSP) {
        if (!fsp) {
            FLOG(ERROR) << "final-state particle empty";
            throw exceptions::Exception("FinalStateParticle empty", "InitialStateParticle::setFinalStateParticles");
        }
        if (!fsp->particleCombinations().empty()) {
            FLOG(ERROR) << "final-state particle already has indices assigned: " << *fsp;
            throw exceptions::Exception("FinalStateParticle already used", "InitialStateParticle::setFinalStateParticles");
        }
        if (fsp->initialStateParticle() != nullptr and fsp->initialStateParticle() != this) {
            FLOG(ERROR) << "final-state particle already has ISP assigned: " << *fsp;
            throw exceptions::Exception("FinalStateParticle already has ISP set", "InitialStateParticle::setFinalStateParticles");
        }
    }

    FinalStateParticles_.reserve(FSP.size());

    // set indices by order in vector
    for (auto& fsp : FSP) {
        fsp->SymmetrizationIndices_.push_back(ParticleCombinationCache_.fsp(FinalStateParticles_.size()));
        fsp->setInitialStateParticle(this);
        FinalStateParticles_.push_back(fsp);
    }
}

//-------------------------
std::array<double, 2> InitialStateParticle::getMassRange(const std::shared_ptr<ParticleCombination>& pc) const
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
ComplexParameterVector InitialStateParticle::freeAmplitudes() const
{
    ComplexParameterVector V = DecayingParticle::freeAmplitudes();

    // if only one free amplitude, treat as fixed and return empty vector
    if (V.size() == 1)
        V.clear();

    return V;
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

    for (unsigned i = 0; i < DataPartitions_.size(); ++i)
        DataPartitions_[i]->setIndex(i);

    setNumberOfDataPartitions(DataPartitions_.size());
}

//-------------------------
void InitialStateParticle::addDataPoint(const std::vector<FourVector<double> >& fourMomenta)
{
    if (!Prepared_) {
        FLOG(ERROR) << "Cannot add DataPoint to InitialStateParticle. Call InitialStateParticle::prepare() first!";
        throw exceptions::Exception("InitialStateParticle not prepared", "InitialStateParticle::addDataPoint");
    }

    if (DataSet_.empty()) {
        addDataPoint(std::move(DataPoint(fourMomenta)));
        return;
    }

    DataSet_.push_back(DataPoint(DataSet_[0]));

    DataPoint& d = DataSet_.back();

    d.setFinalStateFourMomenta(fourMomenta);

    /// calculate FourMomenta_, MeasuredBreakupMomenta_ and HelicityAngles_
    calculate(d);

    if (!DataSet_.consistent(d))
        throw exceptions::Exception("DataPoint inconsistent", "InitialStateParticle::addDataPoint");
}

//-------------------------
void InitialStateParticle::addDataPoint(DataPoint&& d)
{
    if (!Prepared_) {
        LOG(ERROR) << "Cannot add DataPoint to InitialStateParticle. Call InitialStateParticle::prepare() first!";
        throw exceptions::Exception("InitialStateParticle not prepared", "InitialStateParticle::addDataPoint");
    }

    d.allocateStorage(FourMomenta_, DataAccessors_);

    /// calculate FourMomenta_, MeasuredBreakupMomenta_ and HelicityAngles_
    calculate(d);

    if (!DataSet_.consistent(d))
        throw exceptions::Exception("DataPoint inconsistent", "InitialStateParticle::addDataPoint");

    DataSet_.push_back(d);
}

//-------------------------
void InitialStateParticle::addDataPoint(const DataPoint& d)
{ addDataPoint(std::move(DataPoint(d))); }

//-------------------------
void InitialStateParticle::initializeForMonteCarloGeneration(unsigned n)
{
    // initialize with 0
    std::vector<FourVector<double> > momenta(FinalStateParticles_.size(), FourVector_0);

    // add n (empty) data points
    for (unsigned i = 0; i < n; ++i)
        addDataPoint(momenta);

    // set data partitions (1 for each data point)
    setDataPartitions(createDataPartitionsBlockBySize(DataSet_, 1));

    // do one initial calculation
    /// \todo Only calculate data-independent values
    if (!std::isfinite(sumOfLogsOfSquaredAmplitudes()))
        throw exceptions::NonfiniteResult();
}

//-------------------------
void InitialStateParticle::printDataAccessors(bool printParticleCombinations)
{
    // header
    std::cout << "DataAccessors of " << name() << "\n" << "index \tnSymIndices \taddress  \tname";
    if (printParticleCombinations)
        std::cout << "\t\tparticleCombinations";
    std::cout << std::endl;

    for (auto& d : DataAccessors_) {
        std::cout << d->index() << "  \t" << d->maxSymmetrizationIndex() + 1 << "  \t\t" << d << "  \t(" << typeid(*d).name() << ")  \t";
        if (std::dynamic_pointer_cast<Particle>(d))
            std::cout << std::dynamic_pointer_cast<Particle>(d)->name();
        else if (std::dynamic_pointer_cast<DecayChannel>(d))
            std::cout << *std::dynamic_pointer_cast<DecayChannel>(d);

        if (printParticleCombinations) {
            std::cout << " \t";

            for (auto& pc : d->particleCombinations())
                std::cout << *pc << ":" << d->symmetrizationIndex(pc) << ";  ";
        }

        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//-------------------------
void InitialStateParticle::setSymmetrizationIndexParents()
{
    unsigned size = particleCombinations()[0]->indices().size();

    if (size <= 1)
        throw exceptions::Exception("Fewer than two daughters", "InitialStateParticle::setSymmetrizationIndexParents");

    clearSymmetrizationIndices();

    // get initial state PCs from set and add them
    for (auto& wpc : ParticleCombinationCache_) {

        if (wpc.expired())
            continue;

        auto pc = wpc.lock();

        if (pc->indices().size() < size)
            continue;

        if (!pc->daughters()[0]->parent())
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

    FLOG(ERROR) << "trying to calculate for a DataPartition, which is not stored in this InitialStateParticle.";
    return false;
}

//-------------------------
void InitialStateParticle::calculate(DataPoint& d)
{
    // call calculate on static data accessors
    for (auto& sda : DataAccessors_)
        if (std::dynamic_pointer_cast<StaticDataAccessor>(sda))
            std::dynamic_pointer_cast<StaticDataAccessor>(sda)->calculate(d);
}

//-------------------------
void InitialStateParticle::setNumberOfDataPartitions(unsigned n)
{
    // call on this object
    DataAccessor::setNumberOfDataPartitions(n);
    // call on all other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d.get() != this)
            d->setNumberOfDataPartitions(n);
}

//-------------------------
void InitialStateParticle::updateGlobalCalculationStatuses()
{
    // call on this
    DataAccessor::updateGlobalCalculationStatuses();
    // call on all other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d.get() != this)
            d->updateGlobalCalculationStatuses();
}

//-------------------------
void InitialStateParticle::resetCalculationStatuses(unsigned dataPartitionIndex)
{
    // call on this
    DataAccessor::resetCalculationStatuses(dataPartitionIndex);
    // call on other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d.get() != this)
            d->resetCalculationStatuses(dataPartitionIndex);
}

//-------------------------
void InitialStateParticle::setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex)
{
    // call on this
    DataAccessor::setCachedDataValueFlagsToUnchanged(dataPartitionIndex);
    // call on other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d.get() != this)
            d->setCachedDataValueFlagsToUnchanged(dataPartitionIndex);
}

//-------------------------
void InitialStateParticle::setParameterFlagsToUnchanged()
{
    // call on this
    DataAccessor::setParameterFlagsToUnchanged();
    // call on other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d.get() != this)
            d->setParameterFlagsToUnchanged();
}

//-------------------------
DataAccessorSet InitialStateParticle::dataAccessors()
{
    // call DecayingParticle's function
    DataAccessorSet V = DecayingParticle::dataAccessors();

    // add
    V.emplace(FourMomenta_);
    V.emplace(MeasuredBreakupMomenta_);
    V.emplace(HelicityAngles_);

    return V;
}

}
