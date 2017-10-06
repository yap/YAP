#include "Model.h"

#include "Attributes.h"
#include "BlattWeisskopf.h"
#include "CalculationStatus.h"
#include "CompensatedSum.h"
#include "DataAccessor.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "DataSet.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "DecayTree.h"
#include "FinalStateParticle.h"
#include "FreeAmplitude.h"
#include "FourMomenta.h"
#include "Group.h"
#include "logging.h"
#include "MassAxes.h"
#include "Parameter.h"
#include "RecalculableDataAccessor.h"
#include "SpinAmplitudeCache.h"
#include "VariableStatus.h"

/// \todo Find better place for this
INITIALIZE_EASYLOGGINGPP

#include <future>

namespace yap {

//-------------------------
ModelComponent::ModelComponent(const DecayTreeVector& dtv, double adm) :
    DecayTrees_(dtv),
    Admixture_(std::make_shared<NonnegativeRealParameter>(adm))
{
    if (DecayTrees_.empty())
        throw exceptions::Exception("DecayTreeVector empty", "ModelComponent::ModelComponent");
    
    if (group(DecayTrees_, by_parent<>(), by_m<>()).size() != 1)
        throw exceptions::Exception("DecayTree's don't have common parent and spin projection", "ModelComponent::ModelComponent");
}
    
//-------------------------
Model::Model(std::unique_ptr<SpinAmplitudeCache> SAC) :
    Locked_(false),
    FourMomenta_(std::make_shared<FourMomenta>(*this)),
    HelicityAngles_(*this)
{
    if (!SAC)
        throw exceptions::Exception("SpinAmplitudeCache unset", "Model::Model");
    if (!SAC->empty())
        throw exceptions::Exception("SpinAmplitudeCache not empty", "Model::Model");

    setCoordinateSystem({ThreeVector<double>({1., 0., 0.}),
                         ThreeVector<double>({0., 1., 0.}),
                         ThreeVector<double>({0., 0., 1.})});

    SAC->setModel(*this);
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
const double intensity(const ModelComponent& c, const DataPoint& d)
{
    return c.admixture()->value() * intensity(c.decayTrees(), d);
}

//-------------------------
const double intensity(const std::vector<ModelComponent>& C, const DataPoint& d)
{
    return std::accumulate(C.begin(), C.end(), 0., [&d](double I, const auto& c) {return I + intensity(c, d);});
}

//-------------------------
// hidden helper function,
// resolves C++ problem related to naming of functions and call to std::async below
const double sum_of_logs_of_intensities(const Model& M, DataPartition& D, double ped)
{
    // calculate components
    M.calculate(D);

    // if pedestal is zero
    if (ped == 0)
        return std::accumulate(D.begin(), D.end(), CompensatedSum<double>(0.),
                               [&](auto& l, const auto& d){return l += log(intensity(M, d));});
    // else
    return std::accumulate(D.begin(), D.end(), CompensatedSum<double>(0.),
                           [&](auto& l, const auto& d){return l += (log(intensity(M, d)) - ped);});
}

//-------------------------
const double sum_of_log_intensity(const Model& M, DataPartition& D, double ped)
{
    if (M.components().empty())
        throw exceptions::Exception("Model has no components", "sum_of_log_intensity");

    return sum_of_logs_of_intensities(M, D, ped);
}

//-------------------------
const double sum_of_log_intensity(const Model& M, DataPartitionVector& DP, double ped)
{
    // if DataPartitionVector is empty
    if (DP.empty())
        throw exceptions::Exception("DataPartitionVector is empty", "sum_of_log_intensity");

    // check no partitions are nullptr
    if (std::any_of(DP.begin(), DP.end(), std::logical_not<DataPartitionVector::value_type>()))
        throw exceptions::Exception("DataPartitionVector contains nullptr", "sum_of_log_intensity");

    // if threading is unnecessary
    if (DP.size() == 1)
        return sum_of_log_intensity(M, *DP[0], ped);

    if (M.components().empty())
        throw exceptions::Exception("Model has no components", "sum_of_log_intensity");

    std::vector<std::future<double> > partial_sums;
    partial_sums.reserve(DP.size());

    // create thread for calculation on each partition
    for (auto& P : DP)
        // since std::async copies its arguments, even if they are supposed to be references, we need to use std::ref and std::cref
        partial_sums.push_back(std::async(std::launch::async, sum_of_logs_of_intensities, std::cref(M), std::ref(*P), ped));

    // wait for each partition to finish calculating
    return std::accumulate(partial_sums.begin(), partial_sums.end(), 0., [](auto l, auto& s) {return l + s.get();});
}

//-------------------------
bool Model::consistent() const
{
    bool C = true;

    C &= ParticleCombinationCache_.consistent();
    C &= SpinAmplitudeCache_->consistent();

    for (const auto& da : dataAccessors())
        C &= da->consistent();

    for (auto& p : initialStates())
        C &= p->consistent();

    return C;
}

//-------------------------
void Model::addParticleCombination(const ParticleCombination& pc)
{
    if (locked())
        throw exceptions::Exception("Model is locked and cannot be modified.", "Model::addParticleCombination");

    // if does not trace up to an ISP, halt
    if (!is_from_initial_state_particle_combination(pc, *this))
        return;

    FourMomenta_->addParticleCombination(pc);

    // call recursively on daughters
    for (auto& d : pc.daughters())
        addParticleCombination(*d);
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
        fsp->addParticleCombination(*ParticleCombinationCache_.fsp(FinalStateParticles_.size()));
        fsp->setModel(this);
        FinalStateParticles_.push_back(fsp);
    }
}

//-------------------------
void Model::setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm) const
{
    fourMomenta()->setFinalStateMomenta(d, P, sm);

    // call calculate on all static data accessors in model
    for (const auto& sda : StaticDataAccessors_)
        sda->calculate(d, sm);
 
    HelicityAngles_.clearCache(sm);
}

//-------------------------
void Model::addInitialState(std::shared_ptr<DecayingParticle> p)
{
    if (!p)
        throw exceptions::Exception("Initial-state particle empty", "Model::addInitialState");
    addInitialState(*p);
}

//-------------------------
void Model::addInitialState(DecayingParticle& p)
{
    if (locked())
        throw exceptions::Exception("Model is locked and cannot be modified.", "Model::addInitialState");
    
    if (p.model() != this)
        throw exceptions::Exception("Initial-state particle does not belong to this model", "Model::addInitialState");

    auto ptr = std::static_pointer_cast<DecayingParticle>(p.shared_from_this());
    
    // if already present, do nothing
    if (std::find(InitialStates_.begin(), InitialStates_.end(), ptr) != InitialStates_.end())
        return;
    
    // add to InitialStates_
    InitialStates_.push_back(ptr);
}

//-------------------------
size_t all_fixed(const std::vector<ModelComponent>& C)
{
    return std::all_of(C.begin(), C.end(), [](const auto& c){return c.admixture()->variableStatus() == VariableStatus::fixed;});
}

//-------------------------
void fix_solitary_free_amplitudes(Model& m)
{
    // loop over sets of FreeAmplitudes grouped by parent state
    for (auto& fas : group(free_amplitudes(m), by_parent<>()))
        // if there is only one free amplitude for the decay of a state
        if (fas.size() == 1) {
            if (!*fas.begin())
                throw exceptions::Exception("FreeAmplitude is nullptr", "fix_solitary_free_amplitudes");
            // fix it
            (*fas.begin())->variableStatus() = VariableStatus::fixed;
        }
}

//-------------------------
std::vector<std::shared_ptr<DecayingParticle> > full_final_state_isp(const Model& M)
{
    DecayingParticleVector isps;
    isps.reserve(M.initialStates().size());

    // collect first those particles for whom all admixture factors are fixed
    for (const auto& isp : M.initialStates())
        if (decays_to_full_final_state(*isp) and all_fixed(filter(M.components(), from(isp))))
            isps.push_back(isp);

    // then collect the rest
    for (const auto& isp : M.initialStates())
        if (decays_to_full_final_state(*isp) and !all_fixed(filter(M.components(), from(isp))))
            isps.push_back(isp);

    return isps;
}

//-------------------------
void Model::setCoordinateSystem(const CoordinateSystem<double, 3>& cs)
{
    if (!is_right_handed(cs))
        throw exceptions::Exception("Coordinate system not right-handed", "Model::setCoordinateSystem");

    CoordinateSystem_ = unit(cs);
}

//-------------------------
void Model::lock()
{
    // if already locked, do nothing
    if (Locked_)
        return;

    // for each isp
    for (auto& isp : InitialStates_) {

        // register it
        isp->registerWithModel();

        // add its particle combinations
        for (auto& pc : isp->particleCombinations())
            addParticleCombination(*pc);

        // and create components
        for (auto dtv : group(isp->decayTrees(), by_m<>()))
            Components_.emplace_back(dtv);

    }
    
    // call prune() for all initial states
    for (auto& isp : InitialStates_)
        isp->prune();

    fix_solitary_free_amplitudes(*this);

    // if only one ISP with only one spin projection, fix it's admixture
    if (Components_.size() == 1)
        Components_[0].admixture()->variableStatus() = VariableStatus::fixed;

    // remove expired elements of DataAccessors_
    remove_expired(DataAccessors_);
    remove_expired(StaticDataAccessors_);

    // prune remaining DataAccessor's
    for (auto& D : DataAccessors_)
        D->pruneSymmetrizationIndices();

    // remove data accessors from list that don't need storage
    for (auto it = DataAccessors_.begin(); it != DataAccessors_.end(); ) {
        if (!(*it)->requiresStorage())
            it = DataAccessors_.erase(it);
        else
            ++it;
    }

    // set DataAccessor indices
    int index = -1;
    for (const auto& da : DataAccessors_)
        da->setIndex(++index);

    Locked_ = true;
}

//-------------------------
const MassAxes Model::massAxes(std::vector<std::vector<unsigned> > pcs) const
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
    if (std::any_of(pcs.begin(), pcs.end(), [](const auto& v) {return v.size() != 2;}))
        throw exceptions::Exception("only 2-particle axes supported currently", "Model::massAxes");

    ParticleCombinationVector M;

    for (auto& v : pcs) {

        // check that all indices are in range
        {
            auto it = std::find_if(v.end(), v.end(), std::bind(std::greater_equal<unsigned>(), std::placeholders::_1, n_fsp));
            if (it != v.end())
                throw exceptions::Exception("particle index out of range (" + std::to_string(*it) + " >= " + std::to_string(n_fsp) + ")",
                                            "Model::massAxes");
        }

        // check for duplicates
        {
            auto it = std::find_if(v.begin(), v.end(), [&](unsigned i) {return std::count(v.begin(), v.end(), i) != 1;});
            if (it != v.end())
                throw exceptions::Exception("duplicate index (" + std::to_string(*it) + ")given", "Model::massAxes");
        }

        // get fsp ParticleCombinations
        ParticleCombinationVector pcv;
        pcv.reserve(v.size());
        std::transform(v.begin(), v.end(), std::back_inserter(pcv), [&](unsigned i) {return particleCombinationCache().fsp(i);});
        auto pc = particleCombinationCache().composite(pcv);

        // check that pc isn't already in M
        if (std::any_of(M.begin(), M.end(), std::bind(&equal_by_orderless_content, pc, std::placeholders::_1)))
            throw exceptions::Exception("axis requested twice", "Model::massAxes");

        M.push_back(pc);
    }

    return MassAxes(M);
}

//-------------------------
DataSet Model::createDataSet(size_t n)
{
    if (!locked())
        lock();

    if (!locked())
        throw exceptions::Exception("data sets cannot be generated from an unlocked model.", "Model::createDataSet");

    // create empty data set
    DataSet D(*this);

    D.addEmptyDataPoints(n);

    return D;
}

//-------------------------
FreeAmplitudeSet free_amplitudes(const Model& M)
{
    FreeAmplitudeSet S;
    for (const auto& isp : M.initialStates()) {
        auto s = free_amplitudes(*isp);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
ParticleSet particles(const Model& M)
{
    ParticleSet S;
    for (const auto& isp : M.initialStates()) {
        auto s = particles(*isp);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
void Model::setParameterFlagsToUnchanged()
{
    for (auto& d : RecalculableDataAccessors_)
        d->setParameterFlagsToUnchanged();
}

}
