#include "Model.h"

#include "CalculationStatus.h"
#include "Constants.h"
#include "DataAccessor.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "DataSet.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "MassAxes.h"
#include "MeasuredBreakupMomenta.h"
#include "Parameter.h"
#include "PhaseSpaceUtilities.h"
#include "RecalculableDataAccessor.h"
#include "RequiresHelicityAngles.h"
#include "RequiresMeasuredBreakupMomenta.h"
#include "SpinAmplitudeCache.h"
#include "VariableStatus.h"

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
    SAC->setModel(this);
    SpinAmplitudeCache_ = std::move(SAC);
}

//-------------------------
void Model::calculate(DataPartition& D) const
{
    // update calculation statuses
    for (const auto& rda : RecalculableDataAccessors_)
        rda->updateCalculationStatus(D);

    // call calculate on all RecalculableDataAccessors
    for (const auto& rda : RecalculableDataAccessors_)
        rda->calculate(D);
}

//-------------------------
const double intensity(const InitialStateParticleMap& isp_map, const DataPoint& d)
{
    return std::accumulate(isp_map.begin(), isp_map.end(), 0.,
                           [&](double & a, const InitialStateParticleMap::value_type & isp_b2)
    {return a += isp_b2.second->value() * intensity(isp_b2.first->decayTrees(), d);});
}

//-------------------------
// hidden helper function,
// resolves C++ problem related to naming of functions and call to std::async below
const double sum_of_logs_of_intensities(const Model& M, DataPartition& D)
{
    // calculate components
    M.calculate(D);

    // sum log of intensities over data points in partition
    return std::accumulate(D.begin(), D.end(), 0., [&](double & l, const DataPoint & d) {return l += log(intensity(M.initialStateParticles(), d));});
}

//-------------------------
const double sum_of_log_intensity(const Model& M, DataPartition& D)
{
    if (M.initialStateParticles().empty())
        throw exceptions::Exception("Model has no initialStateParticles", "sum_of_log_intensity");

    return sum_of_logs_of_intensities(M, D);
}

//-------------------------
const double sum_of_log_intensity(const Model& M, DataPartitionVector& DP)
{
    // if DataPartitionVector is empty
    if (DP.empty())
        throw exceptions::Exception("DataPartitionVector is empty", "sum_of_log_intensity");

    // check no partitions are nullptr
    if (std::any_of(DP.begin(), DP.end(), std::logical_not<DataPartitionVector::value_type>()))
        throw exceptions::Exception("DataPartitionVector contains nullptr", "sum_of_log_intensity");

    // if threading is unnecessary
    if (DP.size() == 1)
        return sum_of_log_intensity(M, *DP[0]);

    if (M.initialStateParticles().empty())
        throw exceptions::Exception("Model has no InitialStateParticles", "sum_of_log_intensity");

    std::vector<std::future<double> > partial_sums;
    partial_sums.reserve(DP.size());

    // create thread for calculation on each partition
    for (auto& P : DP)
        // since std::async copies its arguments, even if they are supposed to be references, we need to use std::ref and std::cref
        partial_sums.push_back(std::async(std::launch::async, sum_of_logs_of_intensities, std::cref(M), std::ref(*P)));

    // wait for each partition to finish calculating
    return std::accumulate(partial_sums.begin(), partial_sums.end(), 0.,
    [](double & l, std::future<double>& s) {return l += s.get();});
}

//-------------------------
bool Model::consistent() const
{
    bool C = true;

    C &= ParticleCombinationCache_.consistent();
    C &= SpinAmplitudeCache_->consistent();

    for (const auto& da : dataAccessors())
        C &= da->consistent();

    for (auto& p : initialStateParticles())
        C &= p.first->consistent();

    return C;
}

//-------------------------
void Model::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    if (locked())
        throw exceptions::Exception("Model is locked and cannot be modified.", "Model::addParticleCombination");

    // if does not trace up to an ISP, halt
    if (not is_initial_state_particle_combination(*pc, this))
        return;

    FourMomenta_->addParticleCombination(pc);

    if (pc->daughters().size() == 2) {
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
void Model::setFinalState(const FinalStateParticleVector& FSP)
{
    // check that FinalStateParticles_ is empty
    if (!FinalStateParticles_.empty())
        throw exceptions::Exception("Final-state particles already set", "Model::setFinalState");

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

    FinalStateParticles_.reserve(FSP.size());

    // set indices by order in vector
    for (auto& fsp : FSP) {
        fsp->addParticleCombination(ParticleCombinationCache_.fsp(FinalStateParticles_.size()));
        fsp->setModel(this);
        FinalStateParticles_.push_back(fsp);
    }
}

//-------------------------
const InitialStateParticleMap::value_type& Model::addInitialStateParticle(std::shared_ptr<DecayingParticle> p)
{
    if (locked())
        throw exceptions::Exception("Model is locked and cannot be modified.", "Model::addInitialStateParticle");

    if (!p)
        throw exceptions::Exception("Initial-state particle empty", "Model::addInitialStateParticle");

    if (p->model() != this)
        throw exceptions::Exception("Initial-state particle does not belong to this model", "Model::addInitialStateParticle");

    auto res = InitialStateParticles_.insert(std::make_pair(p, std::make_shared<RealParameter>(1.)));

    if (res.second) { // new element was inserted
        for (auto& pc : p->particleCombinations())
            addParticleCombination(pc);

    } else { // no new element was inserted
        // check if insertion failed
        if (res.first == InitialStateParticles_.end())
            throw exceptions::Exception("Failed to insert initialStateParticle", "Model::addInitialStateParticle");

        DEBUG("initial state particle " << *p << " already in model");
    }

    return *res.first;
}

//-------------------------
std::vector<std::shared_ptr<DecayingParticle> > full_final_state_isp(const Model& M)
{
    std::vector<std::shared_ptr<DecayingParticle> > isps;
    isps.reserve(M.initialStateParticles().size());
    for (const auto& kv : M.initialStateParticles())
        if (kv.second->variableStatus() == VariableStatus::fixed
                and kv.first->finalStateParticles(0).size() == M.finalStateParticles().size())
            isps.push_back(kv.first);
    for (const auto& kv : M.initialStateParticles())
        if (kv.second->variableStatus() != VariableStatus::fixed
                and kv.first->finalStateParticles(0).size() == M.finalStateParticles().size())
            isps.push_back(kv.first);
    return isps;
}

//-------------------------
void Model::setCoordinateSystem(const CoordinateSystem<double, 3>& cs)
{
    if (!isRightHanded(cs))
        throw exceptions::Exception("Coordinate system not right-handed", "Model::setCoordinateSystem");

    CoordinateSystem_ = unit(cs);
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

        // if RecalculableDataAccessor, add to RecalculableDataAccessors_
        if (dynamic_cast<RecalculableDataAccessor*>(da))
            RecalculableDataAccessors_.insert(static_cast<RecalculableDataAccessor*>(da));
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

    for (auto& kv : InitialStateParticles_) {
        // prune initial state particles
        kv.first->pruneParticleCombinations();

        // fix amplitudes when they are for the only possible decay chain
        kv.first->fixSolitaryFreeAmplitudes();
    }


    // fix indices:
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
    // if no axes requested, build default:
    if (pcs.empty()) {

        if (finalStateParticles().size() > 4)
            throw exceptions::Exception("Currently only supports final states of 4 or fewer particles for default axes", "Model::massAxes");

        // builds vector down first off diagonal, then second off-diagonal, etc
        for (unsigned i = 1; i < finalStateParticles().size() - 1; ++i)
            for (unsigned j = 0 ; i + j < finalStateParticles().size(); ++j)
                pcs.push_back({j, j + i});
    }

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
        if (std::any_of(v.begin(), v.end(), std::bind(std::greater_equal<unsigned>(), std::placeholders::_1, n_fsp)))
            throw exceptions::Exception("particle index out of range", "Model::massAxes");

        // check for duplicates
        if (std::any_of(v.begin(), v.end(), [&](unsigned i) {return std::count(v.begin(), v.end(), i) != 1;}))
        throw exceptions::Exception("duplicate index given", "Model::massAxes");

        // get fsp ParticleCombinations
        ParticleCombinationVector pcv;
        pcv.reserve(v.size());
        std::transform(v.begin(), v.end(), std::back_inserter(pcv), [&](unsigned i) {return particleCombinationCache().fsp(i);});
        auto pc = particleCombinationCache().composite(pcv);

        // check that pc isn't already in M
        if (any_of(M, pc, ParticleCombination::equalByOrderlessContent))
            throw exceptions::Exception("axis requested twice", "Model::massAxes");

        M.push_back(pc);
    }

    return MassAxes(M);
}

//-------------------------
bool check_invariant_masses(const MassAxes& axes, const std::vector<double>& squared_masses, const std::vector<FourVector<double> >& fourMomenta)
{
    for (size_t i = 0; i < axes.size(); ++i) {
        auto p = std::accumulate(axes[i]->indices().begin(), axes[i]->indices().end(), FourVector_0, [&](const FourVector<double>& p, unsigned j) {return p + fourMomenta[j];});
        if (fabs(norm(p) - squared_masses[i]) > 5. * std::numeric_limits<double>::epsilon() )
            return false;
    }
    return true;
}

//-------------------------
std::vector<FourVector<double> > Model::calculateFourMomenta(const MassAxes& axes, const std::vector<double>& squared_masses, double initial_mass) const
{
    if (FinalStateParticles_.empty())
        throw exceptions::Exception("Final state unset", "Model::calculateFourMomenta");

    if (axes.size() != squared_masses.size())
        throw exceptions::Exception("Incorrect number of masses provided ("
                                    + std::to_string(squared_masses.size()) + " != " + std::to_string(axes.size()) + ")",
                                    "Model::calculateFourMomenta");

    // check none are negative
    if (std::any_of(squared_masses.begin(), squared_masses.end(), std::bind(std::less<double>(), std::placeholders::_1, 0)))
        throw exceptions::Exception("negative squared mass given", "Model::calculateFourMomenta");

    unsigned n_fsp = finalStateParticles().size();

    /// \todo: check sign determination on x component for particles 5 and higher
    if (n_fsp > 4)
        throw exceptions::Exception("not yet supporting 5 or more particles", "Model::calculateFourMomenta");

    if (initial_mass < 0) {
        auto ISPs = full_final_state_isp(*this);
        if (ISPs.empty())
            throw exceptions::Exception("No proper initial-state particle found", "Model::calculateFourMomenta");
        initial_mass = ISPs[0]->mass()->value();
    }

    // matrix of four-vector products
    std::vector<std::vector<double> > pp(n_fsp, std::vector<double>(n_fsp, -1));

    // set diagonal elements to squared masses, and store sum
    double m2_sum_1 = 0;
    double m_sum_1 = initial_mass;
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
                pp[i][j] = pow(initial_mass, 2) + (n_fsp - 2.) * m2_sum_1 - m2_sum_2;
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
            // equal to: if m^2_ij < (m_i + m_j)^2
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
                if (n_fsp == 3 and P[0][3] == 0 and !check_invariant_masses(axes, squared_masses, P))
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
    for (auto& d : RecalculableDataAccessors_)
        d->setParameterFlagsToUnchanged();
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

    for (const auto& d : DataAccessors_) {
        std::cout << d->index() << "  \t" << d->nSymmetrizationIndices() << "  \t\t" << d << "  \t(" << typeid(*d).name() << ")  \t";

        if (printParticleCombinations) {
            std::cout << " \t";

            for (const auto& pc_i : d->symmetrizationIndices())
                std::cout << *pc_i.first << ":" << pc_i.second << ";  ";
        }

        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//-------------------------
void Model::printFlags(const StatusManager& sm) const
{
    for (const auto& d : DataAccessors_) {
        std::cout << std::endl;

        for (auto& c : d->cachedDataValues()) {
            std::cout << "  CachedDataValue " << c << ": ";
            for (unsigned i = 0; i < d->nSymmetrizationIndices(); ++i)
                std::cout << sm.status(*c, i) << "; ";
            std::cout << "\n";
        }
    }

    std::cout << std::endl;
}

}
