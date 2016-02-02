#include "InitialStateParticle.h"

#include "Constants.h"
#include "DataPartition.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "LorentzTransformation.h"
#include "MassAxes.h"
#include "ParticleCombinationCache.h"

#include <assert.h>
#include <future>
#include <memory>
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

    if (!isRightHanded(CoordinateSystem_)) {
        FLOG(ERROR) << "Coordinate system is not right handed.";
        C &= false;
    }

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

    // prepare FourMomenta. Needs FinalStateParticles_
    FourMomenta_->prepare();

    // remove expired elements of DataAccessors_
    removeExpired(DataAccessors_);

    for (auto& D : DataAccessors_)
        D->pruneSymmetrizationIndices();

    FLOG(INFO) << ParticleCombinationCache_;

    for (auto& D : DataAccessors_) {
        std::cout << std::endl;
        D->printParticleCombinations();
    }

    // // check
    // for (auto& wpc : ParticleCombinationCache_) {
    //     if (!wpc.lock())
    //         continue;
    //     auto pc = wpc.lock();
    //     if (!pc->consistent() or pc->origin().indices().size() != finalStateParticles().size())
    //         FLOG(ERROR) << "Cannot prepare InitialStateParticle, particleCombinationCache is not consistent.";
    //     throw exceptions::Exception("ParticleCombination inconsistent", "InitialStateParticle::prepare");
    // }


    // add this (commented out because ISP has no need for data access at moment)
    // DataAccessors_.insert(shared_from_this());
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
void InitialStateParticle::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    // if pc is for this initial state particle, call DecayingParticle::addParticleCombination
    if (pc->indices().size() == FinalStateParticles_.size())
        DecayingParticle::addParticleCombination(pc);

    // if final state particle, halt
    if (pc->isFinalStateParticle())
        return;

    // find top-most parent
    auto p = pc;
    while (p->parent())
        p = p->parent();
    // if does not trace up to ISP, halt
    if (p->indices().size() != FinalStateParticles_.size())
        return;

    FourMomenta_->addParticleCombination(pc);
    HelicityAngles_->addParticleCombination(pc);
    MeasuredBreakupMomenta_->addParticleCombination(pc);

    // call recursively on daughters
    for (auto& d : pc->daughters())
        addParticleCombination(d);
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
        fsp->addParticleCombination(ParticleCombinationCache_.fsp(FinalStateParticles_.size()));
        fsp->setInitialStateParticle(this);
        FinalStateParticles_.push_back(fsp);
    }
}

//-------------------------
void InitialStateParticle::setCoordinateSystem(const CoordinateSystem<double, 3>& cs)
{
    if (!isRightHanded(cs))
        throw exceptions::Exception("Coordinate system not right-handed", "InitialStateParticle::setCoordinateSystem");

    CoordinateSystem_ = unit(cs);
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
    if (!Prepared_)
        throw exceptions::Exception("InitialStateParticle not yet prepared", "InitialStateParticle::addDataPoint");

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
    if (!Prepared_)
        throw exceptions::Exception("ISP not yet prepared", "InitialStateParticle::initializeForMonteCarloGeneration");

    if (!DataSet_.empty())
        throw exceptions::Exception("DataSet isn't empty", "InitialStateParticle::initializeForMonteCarloGeneration");

    // initialize with 0
    std::vector<FourVector<double> > momenta(FinalStateParticles_.size(), FourVector_0);

    // create data point
    auto d = DataPoint(momenta);
    // and allocate space
    d.allocateStorage(FourMomenta_, DataAccessors_);

    // add n (empty) data points
    for (unsigned i = 0; i < n; ++i)
        DataSet_.emplace_back(DataPoint(d));

    // set data partitions (1 for each data point)
    setDataPartitions(createDataPartitionsBlocksBySize(DataSet_, 1));

    // do one initial calculation
    /// \todo Only calculate data-independent values
    // if (!std::isfinite(sumOfLogsOfSquaredAmplitudes()))
    //     throw exceptions::NonfiniteResult();
}

//-------------------------
const MassAxes InitialStateParticle::getMassAxes(std::vector<std::vector<ParticleIndex> > pcs)
{
    unsigned n_fsp = finalStateParticles().size();
    unsigned n_axes = 3 * n_fsp - 7;

    // check that number of requested axes == n_axes
    if (pcs.size() != n_axes) {
        if (pcs.size() < n_axes)
            throw exceptions::Exception("too few axes requested ( " + std::to_string(pcs.size()) + " < " + std::to_string(n_axes) + " )",
                                        "InitialStateParticle::getMassAxes");
        else
            throw exceptions::Exception("too many axes requested ( " + std::to_string(pcs.size()) + " > " + std::to_string(n_axes) + " )",
                                        "InitialStateParticle::getMassAxes");
    }

    // for the moment, we only support 2-particle axes
    // check that all axes are 2 -particle
    if (std::any_of(pcs.begin(), pcs.end(), [](const std::vector<ParticleIndex>& v) {return v.size() != 2;}))
    throw exceptions::Exception("only 2-particle axes supported currently", "InitialStateParticle::getMassAxes");

    ParticleCombinationVector M;

    for (auto& v : pcs) {

        // check that all indices are in range
        if (std::any_of(v.begin(), v.end(), [&](const ParticleIndex & i) {return i >= n_fsp;}))
        throw exceptions::Exception("particle index out of range", "InitialStateParticle::getMassAxes");

        // sort v
        sort(v.begin(), v.end());
        // check for duplicates
        if (std::adjacent_find(v.begin(), v.end()) != v.end())
            throw exceptions::Exception("duplicate index given", "InitialStateParticle::getMassAxes");

        // get ParticleCombination
        auto pc0 = particleCombinationCache().fsp(v[0]);
        auto pc1 = particleCombinationCache().fsp(v[1]);
        auto pc = particleCombinationCache().composite({pc0, pc1});

        // check that pc isn't already in M
        for (const auto& m : M)
            if (ParticleCombination::equivByOrderlessContent(m, pc))
                throw exceptions::Exception("axes requested twice: " + indices_string(*m) + " == " + indices_string(*pc), "InitialStateParticle::getMassAxes");

        M.push_back(pc);
    }

    return MassAxes(M);
}

//-------------------------
bool InitialStateParticle::setMasses(DataPoint& d, const MassAxes& axes, const std::vector<double>& masses)
{
    std::vector<double> squared_masses(masses.size(), -1);
    std::transform(masses.begin(), masses.end(), squared_masses.begin(), [](double m) {return m * m;});
    return setSquaredMasses(d, axes, squared_masses);
}

//-------------------------
bool InitialStateParticle::setSquaredMasses(DataPoint& d, const MassAxes& axes, const std::vector<double>& squared_masses)
{
    if (axes.size() != squared_masses.size())
        throw exceptions::Exception("Incorrect number of masses provided ("
                                    + std::to_string(squared_masses.size()) + " != " + std::to_string(axes.size()) + ")",
                                    "InitialStateParticle::setSquaredMasses");

    // check none are negative
    if (std::any_of(squared_masses.begin(), squared_masses.end(), [](double m) {return m < 0;}))
    throw exceptions::Exception("negative squared mass given", "InitialStateParticle::setSquaredMasses");

    // reset all masses to -1
    FourMomenta_->resetMasses(d);

    unsigned n_fsp = finalStateParticles().size();

    // set two-particle invariant masses for those provided
    std::vector<std::vector<double> > m2_2(n_fsp - 1, std::vector<double>(n_fsp, -1));
    for (size_t i = 0; i < axes.size(); ++i)
        m2_2[axes[i]->indices()[0]][axes[i]->indices()[1]] = squared_masses[i];

    /// \todo: check sign determination on x component for particles 5 and higher
    if (n_fsp > 4)
        throw exceptions::Exception("not yet supporting 5 or more particles", "InitialStateParticle::setSquaredMasses");

    //////////////////////////////////////////////////
    // find unset two-particle invariant mass
    // and calculate fsp squared masses
    unsigned I = 0;
    unsigned J = 1;
    std::vector<double> m2_1(n_fsp, -1);
    for (unsigned i = 0; i < n_fsp; ++i) {
        m2_1[i] = pow(finalStateParticles()[i]->mass()->value(), 2);
        if (m2_2[I][J] >= 0)
            for (unsigned j = i + 1; j < m2_2[i].size(); ++j)
                if (m2_2[i][j] < 0) {
                    I = i;
                    J = j;
                }
    }

    // calculate unset mass: m^2_missing = M^2_isp + (n_fsp - 2) * sum_fsp (m^2) - sum_given(m^2)
    m2_2[I][J] = pow(mass()->value(), 2)
                 + (n_fsp - 2) * std::accumulate(m2_1.begin(), m2_1.end(), 0)
                 - std::accumulate(squared_masses.begin(), squared_masses.end(), 0);

    //////////////////////////////////////////////////
    // calculate all four momenta in m_01 rest frame:
    std::vector<FourVector<double> > P;
    P.reserve(n_fsp);

    // get coordinate system
    auto C = coordinateSystem();

    double m_01 = sqrt(m2_2[0][1]);

    // define p_0 in direction of z axis
    double E0 = (m2_2[0][1] - m2_1[1] + m2_1[0]) / 2. / m_01;
    if (E0 < finalStateParticles()[0]->mass()->value())
        return false;

    double P0 = sqrt(E0 * E0 - m2_1[0]);
    P.push_back(FourVector<double>(E0, P0 * C[2]));

    // define p_1 in direction opposite p_0, with same 3-momentum as p_0
    double E1 = (m2_2[0][1] - m2_1[0] + m2_1[1]) / 2. / m_01;
    if (E1 < finalStateParticles()[1]->mass()->value())
        return false;

    P.push_back(FourVector<double>(E1, -P0 * C[2]));

    // define p_2 to lie in the y-z plane
    double p0p2 = (m2_2[0][2] - m2_1[0] - m2_1[2]) / 2.;
    double p1p2 = (m2_2[1][2] - m2_1[1] - m2_1[2]) / 2.;

    double E2 = (p0p2 + p1p2) / m_01;
    if (E2 < finalStateParticles()[2]->mass()->value())
        return false;

    double P2z = (E0 * p1p2 - E1 * p0p2) / 2. / m_01 / P0;
    double P2y = sqrt(E2 * E2 - m2_1[2] - P2z * P2z);
    P.push_back(FourVector<double>(E2, P2y * C[1] + P2z * C[2]));

    if (n_fsp == 4) {
        // define p_3 to be in positive-x hemisphere
        double p0p3 = (m2_2[0][3] - m2_1[0] - m2_1[3]) / 2.;
        double p1p3 = (m2_2[1][3] - m2_1[1] - m2_1[3]) / 2.;
        double p2p3 = (m2_2[2][3] - m2_1[2] - m2_1[3]) / 2.;

        double E3 = (p0p3 + p1p3) / m_01;
        if (E3 < finalStateParticles()[3]->mass()->value())
            return false;

        double P3z = (E0 * p1p3 - E1 * p0p3) / 2. / m_01 / P0;
        double P3y = (E2 * E3 - p2p3 - P2z * P3z) / P2y;
        double P3x = sqrt(E3 * E3 - m2_1[3] - P3z * P3z - P3y * P3y);
        P.push_back(FourVector<double>(E2, P3x * C[0] + P3y * C[1] + P3z * C[2]));
    }

    // boost:
    // auto p_isp = std::accumulate(P.begin(), P.end(), FourVector_0);
    // auto b = lorentzTransformation<double>(-
    // for (auto& p : P)
    //     p = b * p;

    // hand to data point
    d.setFinalStateFourMomenta(P);
    calculate(d);

    FLOG(INFO) << "out";
    return true;
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
        if (dynamic_cast<Particle*>(d))
            std::cout << dynamic_cast<Particle*>(d)->name();
        else if (dynamic_cast<DecayChannel*>(d))
            std::cout << *dynamic_cast<DecayChannel*>(d);

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
        if (dynamic_cast<StaticDataAccessor*>(sda))
            dynamic_cast<StaticDataAccessor*>(sda)->calculate(d);
}

//-------------------------
void InitialStateParticle::setNumberOfDataPartitions(unsigned n)
{
    // call on this object
    DataAccessor::setNumberOfDataPartitions(n);
    // call on all other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d != this)
            d->setNumberOfDataPartitions(n);
}

//-------------------------
void InitialStateParticle::updateGlobalCalculationStatuses()
{
    // call on this
    DataAccessor::updateGlobalCalculationStatuses();
    // call on all other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d != this)
            d->updateGlobalCalculationStatuses();
}

//-------------------------
void InitialStateParticle::resetCalculationStatuses(unsigned dataPartitionIndex)
{
    // call on this
    DataAccessor::resetCalculationStatuses(dataPartitionIndex);
    // call on other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d != this)
            d->resetCalculationStatuses(dataPartitionIndex);
}

//-------------------------
void InitialStateParticle::setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex)
{
    // call on this
    DataAccessor::setCachedDataValueFlagsToUnchanged(dataPartitionIndex);
    // call on other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d != this)
            d->setCachedDataValueFlagsToUnchanged(dataPartitionIndex);
}

//-------------------------
void InitialStateParticle::setParameterFlagsToUnchanged()
{
    // call on this
    DataAccessor::setParameterFlagsToUnchanged();
    // call on other DataAccessor's (does nothing to StaticDataAccessor's)
    for (auto& d : DataAccessors_)
        if (d != this)
            d->setParameterFlagsToUnchanged();
}

}
