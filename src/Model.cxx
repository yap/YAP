#include "Model.h"

#include "CalculationStatus.h"
#include "Constants.h"
#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "logging.h"
#include "MassAxes.h"
#include "MeasuredBreakupMomenta.h"
#include "PhaseSpaceUtilities.h"
#include "RequiresHelicityAngles.h"
#include "RequiresMeasuredBreakupMomenta.h"
#include "SpinAmplitudeCache.h"

/// \todo Find better place for this
INITIALIZE_EASYLOGGINGPP

#include <future>

namespace yap {

//-------------------------
Model::Model(std::unique_ptr<SpinAmplitudeCache> SAC) :
    Locked_(false),
    CoordinateSystem_(ThreeAxes),
    FourMomenta_(std::make_shared<FourMomenta>(this))
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
std::complex<double> Model::amplitude(DataPoint& d, int two_m, StatusManager& sm) const
{
    if (!InitialStateParticle_)
        throw exceptions::Exception("Initial state unset", "Model::amplitude");

    std::complex<double> a = Complex_0;

    // sum up ISP's amplitude over each particle combination
    for (auto& kv : InitialStateParticle_->symmetrizationIndices()) {
        FDEBUG("calculating for two_m = " << two_m << " and pc = " << *kv.first);
        a += InitialStateParticle_->amplitude(d, kv.first, two_m, sm);
    }

    return a;
}

//-------------------------
std::complex<double> Model::amplitude(DataPoint& d, StatusManager& sm) const
{
    if (!InitialStateParticle_)
        throw exceptions::Exception("Initial state unset", "Model::amplitude");

    if (d.model() != this)
        throw exceptions::Exception("DataPoint is not associated with this model.", "Model::amplitude");

    std::complex<double> a = Complex_0;

    for (int two_m = -InitialStateParticle_->quantumNumbers().twoJ(); two_m <= (int)InitialStateParticle_->quantumNumbers().twoJ(); two_m += 2) {
        for (auto& kv : InitialStateParticle_->symmetrizationIndices()) {
            FDEBUG("calculating for two_m = " << two_m << " and pc = " << *kv.first);
            a += InitialStateParticle_->amplitude(d, kv.first, two_m, sm);
        }
    }

    return a;
}

//-------------------------
double Model::partialSumOfLogsOfSquaredAmplitudes(DataPartitionBase* D, const StatusManager& global) const
{
    double L = 0;

    // loop over data points in partition
    for (DataIterator d = D->begin(); d != D->end(); ++d) {
        DEBUG("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ");
        D->copyCalculationStatuses(global);
        L += log(norm(amplitude(*d, *D)));
    }

    return L;
}

//-------------------------
double Model::sumOfLogsOfSquaredAmplitudes(DataSet& DS) const
{
    // update global calculation statuses (managed by data set)
    DS.globalStatusManager().updateCalculationStatuses(dataAccessors());

    auto log_L = partialSumOfLogsOfSquaredAmplitudes(&DS, DS.globalStatusManager());

    // Set all variable statuses to Unchanged
    DS.globalStatusManager().setAll(VariableStatus::unchanged);

    return log_L;
}


//-------------------------
double Model::sumOfLogsOfSquaredAmplitudes(DataSet& DS, DataPartitionVector& DP) const
{
    // if DataPartitionVector is empty, run over whole data set
    if (DP.empty())
        throw exceptions::Exception("DataPartitionVector is empty", "Model::sumOfLogsOfSquaredAmplitudes");

    // update global calculation statuses (managed by data set)
    DS.globalStatusManager().updateCalculationStatuses(dataAccessors());

    double log_L = 0;

    if (DP.size() == 1) {
        // if threading is unnecessary
        log_L = partialSumOfLogsOfSquaredAmplitudes(DP[0].get(), DS.globalStatusManager());

    } else {

        std::vector<std::future<double> > partial_sums;

        // create thread for calculation on each partition
        for (auto& P : DP)
            // since std::async copies its arguments, even if they are supposed to be references, we need to use std::cref
            partial_sums.push_back(std::async(std::launch::async, &Model::partialSumOfLogsOfSquaredAmplitudes, this, P.get(), std::cref(DS.globalStatusManager())));

        // wait for each partition to finish calculating
        for (auto& s : partial_sums)
            log_L  += s.get();
    }

    // Set all variable statuses to unchanged
    DS.globalStatusManager().setAll(VariableStatus::unchanged);

    // Set all calculation statuses to calculated
    DS.globalStatusManager().setAll(CalculationStatus::calculated);

    return log_L;
}

//-------------------------
bool Model::consistent() const
{
    bool C = true;

    C &= ParticleCombinationCache_.consistent();
    C &= SpinAmplitudeCache_->consistent();

    for (const auto& da : dataAccessors())
        C &= da->consistent();

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
        if (HelicityAngles_)
            HelicityAngles_->addParticleCombination(pc);
        if (MeasuredBreakupMomenta_)
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
std::array<double, 2> Model::massRange(const std::shared_ptr<ParticleCombination>& pc) const
{
    if (!InitialStateParticle_)
        throw exceptions::Exception("Initial state not set", "Model::massRange");

    if (FinalStateParticles_.empty())
        throw exceptions::Exception("Final state not set", "Model::massRange");

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
void Model::addDataAccessor(DataAccessorSet::value_type da)
{
    if (locked())
        throw exceptions::Exception("Model is locked and cannot be modified.", "Model::addDataAccessor");

    // check if already in DataAccessors_
    if (DataAccessors_.find(da) != DataAccessors_.end())
        // do nothing
        return;

    if (da->model() != this)
        throw exceptions::Exception("DataAccessor's Model is not this", "Model::addDataAccessor");

    // add DataAccessor
    if (DataAccessors_.insert(da).second) {
        // if insertion was successful
        // set its index
        da->setIndex(DataAccessors_.size() - 1);

        // if HelicityAngles is empty and DataAccessor requires
        // HelicityAngles, create HelicityAngles. Calling before
        // adding to StaticDataAccessorVector insures that
        // HelicityAngles are called before any newly created
        // StaticDataAccessors
        if (!HelicityAngles_ and dynamic_cast<RequiresHelicityAngles*>(da))
            HelicityAngles_ = std::make_shared<HelicityAngles>(this);

        // if MeasuredBreakupMomenta is empty and DataAccessor required it, create it
        if (!MeasuredBreakupMomenta_ and dynamic_cast<RequiresMeasuredBreakupMomenta*>(da)
                and dynamic_cast<RequiresMeasuredBreakupMomenta*>(da)->requiresMeasuredBreakupMomenta())
            MeasuredBreakupMomenta_ = std::make_shared<MeasuredBreakupMomenta>(this);

        // if StaticDataAccessor, add to StaticDataAccessors_
        if (dynamic_cast<StaticDataAccessor*>(da))
            StaticDataAccessors_.push_back(static_cast<StaticDataAccessor*>(da));
    }
}

//-------------------------
void Model::prepareDataAccessors()
{
    if (locked())
        throw exceptions::Exception("Model is locked and cannot be modified.", "Model::prepareDataAcessors");

    // remove expired elements of DataAccessors_
    removeExpired(DataAccessors_);
    removeExpiredStatic(StaticDataAccessors_);

    // prune remaining DataAccessor's
    for (auto& D : DataAccessors_)
        D->pruneSymmetrizationIndices();

    // fix amplitudes when they are for the only possible decay chain
    for (auto& D : DataAccessors_) {
        if (!dynamic_cast<DecayingParticle*>(D))
            continue;
        // if a decaying particle
        auto C = dynamic_cast<DecayingParticle*>(D)->channels();
        if (C.size() != 1)
            continue;
        // if decaying particle has only one decay channel
        auto A = C[0]->freeAmplitudes().size();
        if (A.size() == 1)
            // if only one free amplitude in decay channel
            A[0]->setVariableStatus(VariableStatus::fixed);
    }

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
                if (da->index() > (int)index)
                    da->setIndex(da->index() - 1);
                used.insert(da->index());
            }
        }

        // if index is now used, increment it by 1
        if (used.find(index) != used.end())
            index += 1;

    }

#ifndef ELPP_DISABLE_DEBUG_LOGS
    std::cout << "StaticDataAccessors:\n";
    for (auto& D : StaticDataAccessors_) {
        std::cout << std::endl;
        D->printParticleCombinations();
    }
    std::cout << "DataAccessors:\n";
    for (auto& D : DataAccessors_) {
        std::cout << std::endl;
        D->printParticleCombinations();
    }
#endif

    lock();
}

//-------------------------
const MassAxes Model::massAxes(std::vector<std::vector<unsigned> > pcs)
{
    unsigned n_fsp = finalStateParticles().size();
    unsigned n_axes = 3 * n_fsp - 7;

    // check that number of requested axes == n_axes
    if (pcs.size() != n_axes) {
        if (pcs.size() < n_axes)
            throw exceptions::Exception("too few axes requested ( " + std::to_string(pcs.size()) + " < " + std::to_string(n_axes) + " )",
                                        "Model::massAxes");
        else
            throw exceptions::Exception("too many axes requested ( " + std::to_string(pcs.size()) + " > " + std::to_string(n_axes) + " )",
                                        "Model::massAxes");
    }

    // for the moment, we only support 2-particle axes
    // check that all axes are 2 -particle
    if (std::any_of(pcs.begin(), pcs.end(), [](const std::vector<unsigned>& v) {return v.size() != 2;}))
    throw exceptions::Exception("only 2-particle axes supported currently", "Model::massAxes");

    ParticleCombinationVector M;

    for (auto& v : pcs) {

        // check that all indices are in range
        if (std::any_of(v.begin(), v.end(), [&](const unsigned & i) {return i >= n_fsp;}))
        throw exceptions::Exception("particle index out of range", "Model::massAxes");

        // sort v
        sort(v.begin(), v.end());
        // check for duplicates
        if (std::adjacent_find(v.begin(), v.end()) != v.end())
            throw exceptions::Exception("duplicate index given", "Model::massAxes");

        // get ParticleCombination
        auto pc0 = particleCombinationCache().fsp(v[0]);
        auto pc1 = particleCombinationCache().fsp(v[1]);
        auto pc = particleCombinationCache().composite({pc0, pc1});

        // check that pc isn't already in M
        for (const auto& m : M)
            if (ParticleCombination::equivByOrderlessContent(m, pc))
                throw exceptions::Exception("axes requested twice: " + indices_string(*m) + " == " + indices_string(*pc), "Model::massAxes");

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
    double m_sum_1 = InitialStateParticle_->mass()->value();
    for (size_t i = 0; i < n_fsp; ++i) {
        pp[i][i] = pow(finalStateParticles()[i]->mass()->value(), 2);
        m2_sum_1 += pp[i][i];
        m_sum_1 -= finalStateParticles()[i]->mass()->value();
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
            // first check if m_ij > (M - sum_m + m_i + m_j)
            if (pp[i][j] > pow(m_sum_1 + finalStateParticles()[i]->mass()->value() + finalStateParticles()[j]->mass()->value(), 2))
                return std::vector<FourVector<double> >();

            // P_i * P_j = (m^2_ij - m^2_i - m^2_j) / 2
            pp[i][j] = (pp[i][j] - pp[i][i] - pp[j][j]) / 2;
            // equivalent to: if m^2_ij < (m_i + m_j)^2
            if (pp[i][j] < finalStateParticles()[i]->mass()->value() * finalStateParticles()[j]->mass()->value())
                return std::vector<FourVector<double> >();

            pp[j][i] = pp[i][j];
        }

    //////////////////////////////////////////////////
    // calculate all four momenta in m_01 rest frame:

    std::vector<FourVector<double> > P(n_fsp, FourVector_0);

    for (unsigned i = 0; i < n_fsp; ++i) {

        FourVector<double> p;

        if (i < 2) {

            P[i][0] = (pp[0][i] + pp[1][i]) / m_01; // E
            P[i][3] = pow_negative_one(i) * P[i][0] * sqrt(1. - pp[i][i] / pow(P[i][0], 2)); // Z

            if (!std::isfinite(P[i][3]))
                return std::vector<FourVector<double> >();

        } else {

            // Energy E_i = (P_0 * P_i + P_1 * P_i) / (E_1 + E_0)
            P[i][0] = (pp[0][i] + pp[1][i]) / (P[0][0] + P[1][0]);

            // if E^2 < m^2
            if (P[i][0] < sqrt(pp[i][i]))
                return std::vector<FourVector<double> >();

            // if p0 and p1 are at rest in m01 rest frame, Z_i := 0
            // else Z_i = (P_1 * P_i - P_0 * P_i - (E_1 - E_0) * E_i) / 2 / Z_0
            P[i][3] = (P[0][3] == 0) ? 0 : (pp[1][i] - pp[0][i] - (P[1][0] - P[0][0]) * P[i][0]) / 2. / P[0][3];

            if (i < 3) {

                // p2 in y-z plane, enforce P^2 = m^2
                P[i][2] = P[i][0] * sqrt(1 - pp[i][i] / pow(P[i][0], 2) - pow(P[i][3] / P[i][0], 2));

                if (!std::isfinite(P[i][2]))
                    return std::vector<FourVector<double> >();

                // phasespace check for special case: p0 and p1 are at rest in m01 rest frame
                if (n_fsp == 3 and P[0][3] == 0 and !checkInvariantMasses(axes, squared_masses, P))
                    return std::vector<FourVector<double> >();

            } else {

                // if Y_2 == 0, Y_i := 0
                // else Y_i = (P_2 * P_i - P_2 * P_i) / Y_2
                P[i][2] = (P[2][2] == 0) ? 0 : (P[2][0] * P[i][0] - pp[2][i] - P[2][3] * P[i][3]) / P[2][2];

                if (i < 4) {

                    // enforce P^2 = m^2
                    P[i][1] = P[i][0] * sqrt(1 - pp[i][i] / pow(P[i][0], 2) - pow(P[i][2] / P[i][0], 2) - pow(P[i][3] / P[i][0], 2));

                    if (!std::isfinite(P[i][1]))
                        return std::vector<FourVector<double> >();

                } else
                    throw exceptions::Exception("not yet supporting 5 or more particles", "Model::calculateFourMomenta");
            }
        }
    }

    // adjust for user-provided coordinate system
    auto C = coordinateSystem();
    for (auto& p : P)
        p = FourVector<double>(p[0], p[1] * C[0] + p[2] * C[1] + p[3] * C[2]);

    return P;
}

//-------------------------
DataSet Model::createDataSet(size_t n)
{
    if (!locked())
        // prepare DataAccessors (locks model)
        prepareDataAccessors();

    if (!locked())
        throw exceptions::Exception("data sets cannot be generated from an unlocked model.", "Model::createDataSet");

    // create empty data set
    DataSet D(*this);

    D.addEmptyPoints(n);

    return D;
}

//-------------------------
void Model::setParameterFlagsToUnchanged()
{
    for (auto& d : DataAccessors_)
        for (auto& c : d->cachedDataValues())
            for (auto& p : c->parameterDependencies())
                if (p->variableStatus() == VariableStatus::changed)
                    p->setVariableStatus(VariableStatus::unchanged);
}

//-------------------------
void Model::printDataAccessors(bool printParticleCombinations) const
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
void Model::printFlags(const StatusManager& sm) const
{
    for (auto& d : DataAccessors_) {
        std::cout << d->data_accessor_type() << "  ";
        if (dynamic_cast<Particle*>(d))
            std::cout << dynamic_cast<Particle*>(d)->name();
        else if (dynamic_cast<DecayChannel*>(d))
            std::cout << *dynamic_cast<DecayChannel*>(d);

        std::cout << std::endl;

        for (auto& c : d->cachedDataValues()) {
            std::cout << "  CachedDataValue " << c << ": ";
            for (int i = 0; i <= d->maxSymmetrizationIndex(); ++i)
                std::cout << sm.status(*c, i) << "; ";
            std::cout << "\n";

            for (auto& p : c->parameterDependencies())
                std::cout << "    depends on Parameter " << p << ": " << p->variableStatus() << "\n";

            for (auto& p : c->cachedDataValueDependencies()) {
                std::cout << "    depends on CachedDataValue " << p << ": ";
                for (int i = 0; i <= p->owner()->maxSymmetrizationIndex(); ++i)
                    std::cout << sm.status(*p, i) << "; ";
                std::cout << "\n";
            }

            for (auto& d : c->daughterCachedDataValueDependencies())
                std::cout << "    depends on daughterCachedDataValue " << d.CDV << ": " << sm.status(*d.CDV, 0) << "; ...\n";
        }
    }

    std::cout << std::endl;
}

//-------------------------
bool Model::checkInvariantMasses(const MassAxes& axes, const std::vector<double>& squared_masses, const std::vector<FourVector<double> >& fourMomenta) const
{
    for (size_t i = 0; i < axes.size(); ++i) {
        auto p = std::accumulate(axes[i]->indices().begin(), axes[i]->indices().end(), FourVector_0, [&](const FourVector<double>& p, unsigned j) {return p + fourMomenta[j];});
        if (fabs(norm(p) - squared_masses[i]) > 5. * std::numeric_limits<double>::epsilon() )
            return false;
    }
    return true;
}

}
