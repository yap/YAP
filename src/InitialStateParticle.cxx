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

    // DEBUG ("InitialStateParticle::amplitude = " << a);

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


    // create data point
    auto d = DataPoint(FinalStateParticles_.size());
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
const MassAxes InitialStateParticle::getMassAxes(std::vector<std::vector<unsigned> > pcs)
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
    if (std::any_of(pcs.begin(), pcs.end(), [](const std::vector<unsigned>& v) {return v.size() != 2;}))
    throw exceptions::Exception("only 2-particle axes supported currently", "InitialStateParticle::getMassAxes");

    ParticleCombinationVector M;

    for (auto& v : pcs) {

        // check that all indices are in range
        if (std::any_of(v.begin(), v.end(), [&](const unsigned & i) {return i >= n_fsp;}))
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

// //-------------------------
// bool InitialStateParticle::setMasses(DataPoint& d, const MassAxes& axes, const std::vector<double>& masses)
// {
//     std::vector<double> squared_masses(masses.size(), -1);
//     std::transform(masses.begin(), masses.end(), squared_masses.begin(), [](double m) {return m * m;});
//     return setSquaredMasses(d, axes, squared_masses);
// }

//-------------------------
std::vector<FourVector<double> > InitialStateParticle::calculateFourMomenta(const MassAxes& axes, const std::vector<double>& squared_masses) const
{
    if (axes.size() != squared_masses.size())
        throw exceptions::Exception("Incorrect number of masses provided ("
                                    + std::to_string(squared_masses.size()) + " != " + std::to_string(axes.size()) + ")",
                                    "InitialStateParticle::setSquaredMasses");

    // check none are negative
    if (std::any_of(squared_masses.begin(), squared_masses.end(), [](double m) {return m < 0;}))
    throw exceptions::Exception("negative squared mass given", "InitialStateParticle::setSquaredMasses");

    unsigned n_fsp = finalStateParticles().size();

    /// \todo: check sign determination on x component for particles 5 and higher
    if (n_fsp > 4)
        throw exceptions::Exception("not yet supporting 5 or more particles", "InitialStateParticle::setSquaredMasses");

    // matrix of four-vector products
    std::vector<std::vector<double> > pp(n_fsp, std::vector<double>(n_fsp, -1));

    // set diagonal elements to squared masses, and store sum
    double m2_sum_1 = 0;
    for (size_t i = 0; i < n_fsp; ++i) {
        pp[i][i] = pow(finalStateParticles()[i]->mass()->value(), 2);
        m2_sum_1 += pp[i][i];
    }

    // add two-particle invariant masses for those provided
    double m2_sum_2 = 0;
    for (size_t i = 0; i < axes.size(); ++i) {
        if (axes[i]->indices()[0] < axes[i]->indices()[1])
            pp[axes[i]->indices()[0]][axes[i]->indices()[1]] = squared_masses[i];
        else
            pp[axes[i]->indices()[1]][axes[i]->indices()[0]] = squared_masses[i];
        m2_sum_2 += squared_masses[i];
    }

    //////////////////////////////////////////////////
    // find unset two-particle invariant mass
    // and calculate fsp squared masses
    for (unsigned i = 0; i < n_fsp; ++i)
        for (unsigned j = i + 1; j < n_fsp; ++j)
            if (pp[i][j] < 0)
                // calculate unset mass: m^2_missing = M^2_isp + (n_fsp - 2) * sum_fsp (m^2) - sum_given(m^2)
                pp[i][j] = pow(mass()->value(), 2) + (n_fsp - 2.) * m2_sum_1 - m2_sum_2;
    // upper triangular elements are now two-particle squared masses

    double m_01 = sqrt(pp[0][1]);

    // finish calculation of off diagonal elements
    for (unsigned i = 0; i < n_fsp; ++i)
        for (unsigned j = i + 1; j < n_fsp; ++j) {
            // if m^2_ij < (m_i + m_j)^2
            if (pp[i][j] < pow(finalStateParticles()[i]->mass()->value() + finalStateParticles()[j]->mass()->value(), 2))
                return std::vector<FourVector<double> >();
            // P_i * P_j = (m^2_ij - m^2_i - m^2_j) / 2
            pp[i][j] = (pp[i][j] - pp[i][i] - pp[j][j]) / 2;
            pp[j][i] = pp[i][j];
        }

    // if (n_fsp == 3) {
    //     double s = 0;
    //     for (unsigned i = 0; i < n_fsp; ++i) {
    //         double j = (i + 1) % n_fsp;
    //         double k = (i + 2) % n_fsp;
    //         s += pow(pp[j][k], 2) / pp[j][j] / pp[k][k];
    //     }
    //     if (s < 1)
    //         return std::vector<FourVector<double> >();
    // }

    //////////////////////////////////////////////////
    // calculate all four momenta in m_01 rest frame:

    std::vector<FourVector<double> > P(n_fsp, FourVector_0);

    for (unsigned i = 0; i < n_fsp; ++i) {

        const double E = (pp[0][i] + pp[1][i]) / m_01;
        const double E2 = E * E;

        // if E^2 < m^2
        if (E2 < pp[i][i])
            return std::vector<FourVector<double> >();

        if (i < 2) {

            const double Z = sqrt(E2 - pp[i][i]);

            // p0 in z direction, p1 in negative z direction
            P[i] = {E, 0, 0, pow_negative_one(i)* Z};

        } else {

            // Z_i = (E_0 * P_1 * P_i - E_1 * P_0 * P_i) / m_01 / Z_0
            const double Z = (P[0][0] * pp[1][i] - P[1][0] * pp[0][i]) / m_01 / P[0][3];

            if (!std::isfinite(Z))
                return std::vector<FourVector<double> >();

            if (i < 3) {

                // p2 in y-z plane
                const double Y = sqrt(E2 - pp[i][i] - Z * Z);

                if (!std::isfinite(Y))
                    return std::vector<FourVector<double> >();

                P[i] = {E, 0, Y, Z};

            } else {

                // Y_i = (E_2 * E_i - P_2 * P_i - Z_2 * Z_i) / Y_2
                const double Y = (P[2][0] * E - pp[2][i] - P[2][3] * Z) / P[2][2];

                if (!std::isfinite(Y))
                    return std::vector<FourVector<double> >();

                if (i < 4) {

                    const double X = sqrt(E2 - pp[3][3] - Z * Z - Y * Y);

                    if (!std::isfinite(X))
                        return std::vector<FourVector<double> >();

                    P[i] = {E, X, Y, Z};

                } else
                    throw exceptions::Exception("not yet supporting 5 or more particles", "InitialStateParticle::setSquaredMasses");
            }
        }
    }

    // adjust for user-provided coordinate system
    auto C = coordinateSystem();
    for (auto& p : P)
        p = FourVector<double>(p[0], p[1] * C[0] + p[2] * C[1] + p[3] * C[2]);

    return P;

    // // calculate energies in m_01 rest frame:
    // std::vector<double> E(n_fsp, -1);
    // for (unsigned i = 0; i < n_fsp; ++i) {
    //     E[i] = (pp[0][i] + pp[1][i]) / m_01;
    //     if (E[i] < finalStateParticles()[i]->mass()->value())
    //         return std::vector<FourVector<double> >();
    // }

    // // calculate Z components
    // std::vector<double> Z(n_fsp, -1);
    // // define p0 in +z direction
    // Z[0] = sqrt(E[0] * E[0] - pp[0][0]);
    // // define p1 in -z direction
    // Z[1] = -Z[0];
    // // rest are calculated from p0 and p1
    // for (unsigned i = 2; i < n_fsp; ++i)
    //     Z[i] = (E[0] * pp[1][i] - E[1] * pp[0][i]) / m_01 / Z[0];

    // // calculate Y components
    // std::vector<double> Y(n_fsp, 0);
    // // p0 and p1 have 0 y component
    // // p2 is defined in z-y plane
    // Y[2] = sqrt(E[2] * E[2] - pp[2][2] - Z[2] * Z[2]);
    // // rest are calculated from p2
    // for (unsigned i = 3; i < n_fsp; ++i)
    //     Y[i] = (E[2] * E[i] - pp[2][i] - Z[2] * Z[i]) / Y[2];

    // // calculate X components
    // std::vector<double> X(n_fsp, 0);
    // // p0, p1, p2 are defined in z-y plane
    // // p3 defines positive x direction
    // if (3 < n_fsp)
    //     Y[3] = sqrt(E[3] * E[3] - pp[3][3] - Z[3] * Z[3] - Y[3] * Y[3]);

    // std::vector<FourVector<double> > P;
    // P.reserve(n_fsp);

    // // get coordinate system
    // auto C = coordinateSystem();
    // for (unsigned i = 0; i < n_fsp; ++i)
    //     P.push_back(FourVector<double>(E[i], X[i] * C[0] + Y[i] * C[1] + Z[i] * C[2]));

    // boost:
    // auto p_isp = std::accumulate(P.begin(), P.end(), FourVector_0);
    // auto b = lorentzTransformation<double>(-
    // for (auto& p : P)
    //     p = b * p;

    // return P;
}

//-------------------------
void InitialStateParticle::setFinalStateFourMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, unsigned dataPartitionIndex)
{
    if (!d.setFinalStateFourMomenta(P, true))
        // if FSP four momenta are changed
        calculate(d, dataPartitionIndex);
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
void InitialStateParticle::calculate(DataPoint& d, unsigned dataPartitionIndex)
{
    // call calculate on static data accessors
    for (auto& sda : DataAccessors_)
        if (dynamic_cast<StaticDataAccessor*>(sda))
            dynamic_cast<StaticDataAccessor*>(sda)->calculate(d, dataPartitionIndex);
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
