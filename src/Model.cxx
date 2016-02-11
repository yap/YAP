#include "Model.h"

#include "Constants.h"
#include "DataPartition.h"
#include "DataSet.h"
#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "LorentzTransformation.h"
#include "MassAxes.h"
#include "ParticleCombinationCache.h"

#include <assert.h>
#include <future>
#include <memory>

namespace yap {

//-------------------------
Model::Model(std::unique_ptr<SpinAmplitudeCache> SAC) :
    CoordinateSystem_(ThreeAxes),
    FourMomenta_(std::make_shared<FourMomenta>(this)),
    MeasuredBreakupMomenta_(std::make_shared<MeasuredBreakupMomenta>(this)),
    HelicityAngles_(std::make_shared<HelicityAngles>(this))
{
    if (!SAC)
        throw exceptions::Exception("SpinAmplitudeCache unset", "Model::Model");
    if (!SAC->empty())
        throw exceptions::Exception("SpinAmplitudeCache not empty", "Model::Model");
    if (SAC->model())
        throw exceptions::Exception("SpinAmplitudeCache already has owner", "Model::Model");
    SpinAmplitudeCache_ = std::move(SAC);
    SpinAmplitudeCache_->setModel(this);
}

//-------------------------
std::complex<double> Model::amplitude(DataPoint& d, int two_m, unsigned dataPartitionIndex) const
{
    if (!InitialStateParticle_)
        throw exceptions::Exception("Initial state unset", "Model::amplitude");

    std::complex<double> a = Complex_0;

    // sum up ISP's amplitude over each particle combination
    for (auto& kv : InitialStateParticle_->symmetrizationIndices()) {
        FDEBUG("calculating for two_m = " << two_m << " and pc = " << *kv.first);
        a += InitialStateParticle_->amplitude(d, kv.first, two_m, dataPartitionIndex);
    }

    return a;
}

//-------------------------
std::complex<double> Model::amplitude(DataPoint& d, unsigned dataPartitionIndex) const
{
    if (!InitialStateParticle_)
        throw exceptions::Exception("Initial state unset", "Model::amplitude");

    std::complex<double> a = Complex_0;

    for (int two_m = -InitialStateParticle_->quantumNumbers().twoJ(); two_m <= (int)InitialStateParticle_->quantumNumbers().twoJ(); two_m += 2) {
        for (auto& kv : InitialStateParticle_->symmetrizationIndices()) {
            FDEBUG("calculating for two_m = " << two_m << " and pc = " << *kv.first);
            a += InitialStateParticle_->amplitude(d, kv.first, two_m, dataPartitionIndex);
        }
    }

    return a;
}

//-------------------------
double Model::logOfSquaredAmplitude(DataPoint& d, unsigned dataPartitionIndex)
{
    resetCalculationStatuses(dataPartitionIndex);
    return log(norm(amplitude(d, dataPartitionIndex)));
}

//-------------------------
double Model::partialSumOfLogsOfSquaredAmplitudes(DataPartitionBase* D)
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
double Model::sumOfLogsOfSquaredAmplitudes()
{
    // update global caclulationStatus's before looping over partitions
    updateGlobalCalculationStatuses();

    // create thread for calculation on each partition
    std::vector<std::future<double> > L;
    for (auto& D : DataPartitions_)
        L.push_back(std::async(std::launch::async, &Model::partialSumOfLogsOfSquaredAmplitudes, this, D.get()));

    // wait for each partition to finish calculating
    double sum = 0;
    for (auto& l : L)
        sum += l.get();

    setParameterFlagsToUnchanged();

    return sum;
}

//-------------------------
bool Model::consistent() const
{
    bool C = true;

    C &= FourMomenta_->consistent();
    C &= MeasuredBreakupMomenta_->consistent();
    C &= HelicityAngles_->consistent();
    C &= ParticleCombinationCache_.consistent();
    C &= SpinAmplitudeCache_->consistent();

    if (!isRightHanded(CoordinateSystem_)) {
        FLOG(ERROR) << "Coordinate system is not right handed.";
        C &= false;
    }

    /// \todo: is this necessary to check?
    /// @{
    std::vector<std::shared_ptr<FinalStateParticle> > A = FinalStateParticles_;
    sort(A.begin(), A.end());
    std::vector<std::shared_ptr<FinalStateParticle> > B = InitialStateParticle_->finalStateParticles();
    sort(B.begin(), B.end());

    if (A != B) {
        FLOG(ERROR) << "FinalStateParticles_ are not set correctly.";
        C &= false;
    }
    /// @}

    C &= InitialStateParticle_->consistent();

    return C;
}

//-------------------------
void Model::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    // find top-most parent
    auto p = pc;
    while (p->parent())
        p = p->parent();
    // if does not trace up to ISP, halt
    if (p->indices().size() != FinalStateParticles_.size())
        return;

    FourMomenta_->addParticleCombination(pc);

    if (!pc->isFinalStateParticle()) {
        HelicityAngles_->addParticleCombination(pc);
        MeasuredBreakupMomenta_->addParticleCombination(pc);
    }

    // call recursively on daughters
    for (auto& d : pc->daughters())
        addParticleCombination(d);
}

//-------------------------
void Model::setInitialStateParticle(std::shared_ptr<DecayingParticle> isp)
{
    if (InitialStateParticle_)
        throw exceptions::Exception("Initial-state particle already set", "Model::setInitialStateParticle");

    if (!isp)
        throw exceptions::Exception("Initial-state particle empty", "Model::setInitialStateParticle");

    if (isp->model() != this)
        throw exceptions::Exception("Initial-state particle does not belong to this model", "Model::setInitialStateParticle");

    if (isp->finalStateParticles().size() != FinalStateParticles_.size())
        throw exceptions::Exception("Initial-state particle has wrong number of final-state particles ("
                                    + std::to_string(isp->finalStateParticles().size()) + " != " + std::to_string(FinalStateParticles_.size()) + ")",
                                    "Model::setInitialStateParticle");

    InitialStateParticle_ = isp;
}

//-------------------------
void Model::setFinalState(std::initializer_list<std::shared_ptr<FinalStateParticle> > FSP)
{
    FDEBUG("1");

    // check that FinalStateParticles_ is empty
    if (!FinalStateParticles_.empty())
        throw exceptions::Exception("Final-state particles already set", "Model::setFinalState");

    FDEBUG("2");

    // check that none of the FSP's has yet been used
    // and that FinalStateParticles don't have ISP set to another ISP
    for (auto& fsp : FSP) {
        if (!fsp)
            throw exceptions::Exception("FinalStateParticle empty", "Model::setFinalState");
        if (!fsp->particleCombinations().empty())
            throw exceptions::Exception("FinalStateParticle already used", "Model::setFinalState");
        if (fsp->model() != nullptr)
            throw exceptions::Exception("FinalStateParticle already has Model set", "Model::setFinalState");
    }

    FDEBUG("3");

    FinalStateParticles_.reserve(FSP.size());

    // set indices by order in vector
    for (auto& fsp : FSP) {
        fsp->addParticleCombination(ParticleCombinationCache_.fsp(FinalStateParticles_.size()));
        fsp->setModel(this);
        FinalStateParticles_.push_back(fsp);
    }
    FDEBUG("4");
}

//-------------------------
void Model::setCoordinateSystem(const CoordinateSystem<double, 3>& cs)
{
    if (!isRightHanded(cs))
        throw exceptions::Exception("Coordinate system not right-handed", "Model::setCoordinateSystem");

    CoordinateSystem_ = unit(cs);
}

//-------------------------
std::array<double, 2> Model::getMassRange(const std::shared_ptr<ParticleCombination>& pc) const
{
    if (!InitialStateParticle_)
        throw exceptions::Exception("Initial state not set", "Model::getMassRange");

    if (FinalStateParticles_.empty())
        throw exceptions::Exception("Final state not set", "Model::getMassRange");

    std::array<double, 2> m = {0, InitialStateParticle_->mass()->value()};

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
ComplexParameterVector Model::freeAmplitudes() const
{
    return InitialStateParticle_->freeAmplitudes();
}

//-------------------------
std::vector<DataPartitionBase*> Model::dataPartitions()
{
    std::vector<DataPartitionBase*> partitions;
    partitions.reserve(DataPartitions_.size());

    for (auto& d : DataPartitions_)
        partitions.push_back(d.get());

    return partitions;
}

//-------------------------
void Model::setDataPartitions(std::vector<std::unique_ptr<DataPartitionBase> > partitions)
{
    DataPartitions_ = std::move(partitions);

    for (unsigned i = 0; i < DataPartitions_.size(); ++i)
        DataPartitions_[i]->setIndex(i);

    setNumberOfDataPartitions(DataPartitions_.size());
}

//-------------------------
void Model::addDataPoint(const std::vector<FourVector<double> >& fourMomenta)
{
    // if adding first data point
    if (DataSet_.empty()) {
        // prepare data accessors
        prepareDataAccessors();
        // and add first point
    }

    DataSet_.emplace_back(DataAccessors_);
    auto& d = DataSet_.back();

#ifndef ELPP_DISABLE_DEBUG_LOGS
    if (DataSet_.size() == 1)
        for (const auto& da : DataAccessors_)
            FDEBUG("assigned  " << da->data_accessor_type() << " at index " << da->index() << " a vector of size " << d.nElements(da->index()));
#endif

    if (!DataSet_.consistent(d))
        throw exceptions::Exception("produced inconsistent data point", "Model::addDataPoint");

    setFinalStateMomenta(d, fourMomenta);
}

//-------------------------
void Model::addDataAccessor(DataAccessorSet::value_type da)
{
    // check if already in DataAccessors_
    if (DataAccessors_.find(da) != DataAccessors_.end())
        // do nothing
        return;

    if (da->model() != this)
        throw exceptions::Exception("DataAccessor's Model is not this", "Model::addDataAccessor");

    if (DataAccessors_.insert(da).second)
        // if insertion was successful
        da->setIndex(DataAccessors_.size() - 1);
}

//-------------------------
void Model::prepareDataAccessors()
{
    // remove expired elements of DataAccessors_
    removeExpired(DataAccessors_);

    // prune remaining DataAccessor's
    for (auto& D : DataAccessors_)
        D->pruneSymmetrizationIndices();

    // fix indices

    // collect used indices
    std::set<unsigned> used;
    for (const auto& da : DataAccessors_)
        used.insert(da->index());

    // repair
    unsigned index = 0;
    while (index < used.size()) {

        // if index is not used
        if (used.find(index) == used.end()) {
            // clear used
            used.clear();
            // reduce all DataAccessor indices greater than index by 1
            // and rebuild used
            for (auto& da : DataAccessors_) {
                if (da->index() > index)
                    da->setIndex(da->index() - 1);
                used.insert(da->index());
            }
        }

        // if index is now used, increment it by 1
        if (used.find(index) != used.end())
            index += 1;

    }

    for (auto& D : DataAccessors_) {
        std::cout << std::endl;
        D->printParticleCombinations();
    }
}

//-------------------------
void Model::initializeForMonteCarloGeneration(unsigned n)
{
    if (!DataSet_.empty())
        throw exceptions::Exception("DataSet isn't empty", "Model::initializeForMonteCarloGeneration");

    prepareDataAccessors();

    // create data point
    auto d = DataPoint(DataAccessors_);

#ifndef ELPP_DISABLE_DEBUG_LOGS
    for (const auto& da : DataAccessors_)
        FDEBUG("assigned  " << da->data_accessor_type() << " at index " << da->index() << " a vector of size " << d.nElements(da->index()));
#endif

    // add n (empty) data points
    for (unsigned i = 0; i < n; ++i)
        DataSet_.emplace_back(DataPoint(d));

    // set data partitions (1 for each data point)
    setDataPartitions(createDataPartitionsBlocksBySize(DataSet_, 1));
}

//-------------------------
const MassAxes Model::getMassAxes(std::vector<std::vector<unsigned> > pcs)
{
    unsigned n_fsp = finalStateParticles().size();
    unsigned n_axes = 3 * n_fsp - 7;

    // check that number of requested axes == n_axes
    if (pcs.size() != n_axes) {
        if (pcs.size() < n_axes)
            throw exceptions::Exception("too few axes requested ( " + std::to_string(pcs.size()) + " < " + std::to_string(n_axes) + " )",
                                        "Model::getMassAxes");
        else
            throw exceptions::Exception("too many axes requested ( " + std::to_string(pcs.size()) + " > " + std::to_string(n_axes) + " )",
                                        "Model::getMassAxes");
    }

    // for the moment, we only support 2-particle axes
    // check that all axes are 2 -particle
    if (std::any_of(pcs.begin(), pcs.end(), [](const std::vector<unsigned>& v) {return v.size() != 2;}))
    throw exceptions::Exception("only 2-particle axes supported currently", "Model::getMassAxes");

    ParticleCombinationVector M;

    for (auto& v : pcs) {

        // check that all indices are in range
        if (std::any_of(v.begin(), v.end(), [&](const unsigned & i) {return i >= n_fsp;}))
        throw exceptions::Exception("particle index out of range", "Model::getMassAxes");

        // sort v
        sort(v.begin(), v.end());
        // check for duplicates
        if (std::adjacent_find(v.begin(), v.end()) != v.end())
            throw exceptions::Exception("duplicate index given", "Model::getMassAxes");

        // get ParticleCombination
        auto pc0 = particleCombinationCache().fsp(v[0]);
        auto pc1 = particleCombinationCache().fsp(v[1]);
        auto pc = particleCombinationCache().composite({pc0, pc1});

        // check that pc isn't already in M
        for (const auto& m : M)
            if (ParticleCombination::equivByOrderlessContent(m, pc))
                throw exceptions::Exception("axes requested twice: " + indices_string(*m) + " == " + indices_string(*pc), "Model::getMassAxes");

        M.push_back(pc);
    }

    return MassAxes(M);
}

//-------------------------
std::vector<FourVector<double> > Model::calculateFourMomenta(const MassAxes& axes, const std::vector<double>& squared_masses) const
{
    if (!InitialStateParticle_)
        throw exceptions::Exception("Initial state unset", "Model::calculateFourMomenta");

    if (FinalStateParticles_.empty())
        throw exceptions::Exception("Final state unset", "Model::calculateFourMomenta");

    if (axes.size() != squared_masses.size())
        throw exceptions::Exception("Incorrect number of masses provided ("
                                    + std::to_string(squared_masses.size()) + " != " + std::to_string(axes.size()) + ")",
                                    "Model::calculateFourMomenta");

    // check none are negative
    if (std::any_of(squared_masses.begin(), squared_masses.end(), [](double m) {return m < 0;}))
    throw exceptions::Exception("negative squared mass given", "Model::calculateFourMomenta");

    unsigned n_fsp = finalStateParticles().size();

    /// \todo: check sign determination on x component for particles 5 and higher
    if (n_fsp > 4)
        throw exceptions::Exception("not yet supporting 5 or more particles", "Model::calculateFourMomenta");

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
                pp[i][j] = pow(InitialStateParticle_->mass()->value(), 2) + (n_fsp - 2.) * m2_sum_1 - m2_sum_2;
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
                    throw exceptions::Exception("not yet supporting 5 or more particles", "Model::calculateFourMomenta");
            }
        }
    }

    // adjust for user-provided coordinate system
    auto C = coordinateSystem();
    for (auto& p : P)
        p = FourVector<double>(p[0], p[1] * C[0] + p[2] * C[1] + p[3] * C[2]);

    // boost:
    // auto p_isp = std::accumulate(P.begin(), P.end(), FourVector_0);
    // auto b = lorentzTransformation<double>(-p_isp);
    // for (auto& p : P)
    //     p = b * p;

    return P;
}

//-------------------------
void Model::setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, unsigned dataPartitionIndex)
{
    FDEBUG("1");
    FourMomenta_->setFinalStateMomenta(d, P, dataPartitionIndex);
    FDEBUG("2");
    calculate(d, dataPartitionIndex);
    FDEBUG("3");
}

//-------------------------
void Model::printDataAccessors(bool printParticleCombinations)
{
    // header
    std::cout << "DataAccessors of \n"
              << "index \tnSymIndices \taddress  \tname";

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
bool Model::hasDataPartition(DataPartitionBase* d)
{
    for (auto& dp : DataPartitions_)
        if (d == dp.get())
            return true;

    return false;
}

//-------------------------
void Model::calculate(DataPoint& d, unsigned dataPartitionIndex)
{
    // these need to be calculated in order!
    FourMomenta_->calculate(d, dataPartitionIndex);
    MeasuredBreakupMomenta_->calculate(d, dataPartitionIndex);
    HelicityAngles_->calculate(d, dataPartitionIndex);

    // \todo find a nicer solution

    // call calculate on static data accessors
    for (auto& sda : DataAccessors_) {
        if (sda == static_cast<StaticDataAccessor*>(FourMomenta_.get()) or
                sda == static_cast<StaticDataAccessor*>(MeasuredBreakupMomenta_.get()) or
                sda == static_cast<StaticDataAccessor*>(HelicityAngles_.get()))
            continue;

        if (dynamic_cast<StaticDataAccessor*>(sda))
            static_cast<StaticDataAccessor*>(sda)->calculate(d, dataPartitionIndex);
    }
}

//-------------------------
void Model::setNumberOfDataPartitions(unsigned n)
{
    for (auto& d : DataAccessors_)
        d->setNumberOfDataPartitions(n);
}

//-------------------------
void Model::updateGlobalCalculationStatuses()
{
    for (auto& d : DataAccessors_)
        d->updateGlobalCalculationStatuses();
}

//-------------------------
void Model::resetCalculationStatuses(unsigned dataPartitionIndex)
{
    for (auto& d : DataAccessors_)
        d->resetCalculationStatuses(dataPartitionIndex);
}

//-------------------------
void Model::setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex)
{
    for (auto& d : DataAccessors_)
        d->setCachedDataValueFlagsToUnchanged(dataPartitionIndex);
}

//-------------------------
void Model::setParameterFlagsToUnchanged()
{
    for (auto& d : DataAccessors_)
        d->setParameterFlagsToUnchanged();
}

}
