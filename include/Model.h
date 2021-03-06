/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

#ifndef yap_Model_h
#define yap_Model_h

#include "fwd/Model.h"

#include "fwd/DataAccessor.h"
#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/DataSet.h"
#include "fwd/DecayingParticle.h"
#include "fwd/DecayTree.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/FourMomenta.h"
#include "fwd/FourVector.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/MassAxes.h"
#include "fwd/Parameter.h"
#include "fwd/Particle.h"
#include "fwd/RecalculableDataAccessor.h"
#include "fwd/StaticDataAccessor.h"
#include "fwd/StatusManager.h"

#include "CoordinateSystem.h"
#include "Filter.h"
#include "HelicityAngles.h"
#include "ParticleCombinationCache.h"
#include "SpinAmplitudeCache.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

/// \class ModelComponent
/// a term to be summed into the intensity incoherently,
/// the term value is Admixture * coherent sum of DecayTrees
class ModelComponent
{
public:
    /// Constructor
    ModelComponent(const DecayTreeVector& dtv, double adm = 1.);

    /// \return DecayTrees (const)
    const DecayTreeVector& decayTrees() const
    { return DecayTrees_; }

    /// \return Admixture_ (const)
    const std::shared_ptr<NonnegativeRealParameter>& admixture() const
    { return Admixture_; }

    /// \return Admixture_
    std::shared_ptr<NonnegativeRealParameter>& admixture()
    { return Admixture_; }

private:
    /// DecayTrees to be coherently summed
    DecayTreeVector DecayTrees_;

    /// incoherent sum admixture (real) multiplying intensity of DecayTrees
    std::shared_ptr<NonnegativeRealParameter> Admixture_;
};

/// \class Model
/// \brief Class implementing a PWA model
/// \author Johannes Rauch, Daniel Greenwald
class Model
{
public:

    /// Constructor
    /// \param SAC unique_ptr to SpinAmplitudeCache
    Model(std::unique_ptr<SpinAmplitudeCache> SAC);

    /// copy constructor deleted
    /// \todo Implement deep copy.
    Model(const Model&) = delete;

    /// copy assignment operator
    /// \todo Implement deep copy assignment.
    Model& operator=(const Model&) = delete;

    /// move constructor
    /// \todo Implement move construction.
    Model(Model&&) = delete;

    /// move assignment operator
    /// \todo Implement move assignment.
    Model& operator=(Model&&) = delete;

    /// Virtual destructor = default
    virtual ~Model() = default;
    
    /// Calculate model for each data point in the data partition
    /// \param D DataPartition to calculate over
    /// \todo This need not be a member function!
    void calculate(DataPartition& D) const;

    /// Check consistency of object
    virtual bool consistent() const;

    /// \return whether model has been locked
    bool locked() const
    { return Locked_; }

    /// prepare model and mark as locked:
    /// removes expired DataAccessor's, prune's remaining, and assigns them indices;
    /// fixes amplitudes that needn't be free
    void lock();

    /// \name Getters
    /// @{

    /// \return coordinate system (const)
    const CoordinateSystem<double, 3>& coordinateSystem() const
    { return CoordinateSystem_; }

    /// \return FourMomenta accessor
    std::shared_ptr<FourMomenta> fourMomenta()
    { return FourMomenta_; }

    /// \return FourMomenta accessor (const)
    const std::shared_ptr<FourMomenta> fourMomenta() const
    { return FourMomenta_; }

    /// \return HelicityAngles accessor
    const HelicityAngles& helicityAngles() const
    { return HelicityAngles_; }

    /// \return ParticleCombinationCache
    ParticleCombinationCache& particleCombinationCache() const
    { return ParticleCombinationCache_; }

    /// \return SpinAmplitudeCache
    SpinAmplitudeCache* spinAmplitudeCache()
    { return SpinAmplitudeCache_.get(); }

    /// \return SpinAmplitudeCache (const)
    const SpinAmplitudeCache* spinAmplitudeCache() const
    { return SpinAmplitudeCache_.get(); }

    /// \return vector of shared pointers to final state particles
    const FinalStateParticleVector& finalStateParticles() const
    { return FinalStateParticles_; }

    /// \return InitialStates_
    const DecayingParticleVector& initialStates() const
    { return InitialStates_; }

    /// \return Components_
    const std::vector<ModelComponent>& components() const
    { return Components_; }
    
    /// \return set of DataAccessors
    const DataAccessorSet& dataAccessors() const
    { return DataAccessors_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set final-state particle content. The order in which particles
    /// are given dictates the order in which four-momenta must be
    /// given in data points. The FinalStateParticle's have their
    /// Model pointer set to this
    /// \param FSP vector of shared pointers to final-state particles
    void setFinalState(const FinalStateParticleVector& FSP);

    /// Set final-state particle content. The order in which particles
    /// are given dictates the order in which four-momenta must be
    /// given in data points. The FinalStateParticle's have their
    /// Model pointer set to this
    /// \param FSPs shared pointers to final-state particles
    template <typename ... Types>
    void setFinalState(Types ... FSPs)
    { FinalStateParticleVector V{FSPs...}; setFinalState(V); }

    /// add an initial state
    void addInitialState(std::shared_ptr<DecayingParticle> p);

    /// add an initial state
    void addInitialState(DecayingParticle& p);

    /// set four momenta of data point
    /// \param P vector of FourVectors of final-state momenta
    /// \param sm StatusManager to update
    void setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm) const;

    /// set coordinate system
    void setCoordinateSystem(const CoordinateSystem<double, 3>& cs);

    /// @}

    /// \name Monte Carlo Generation
    /// @{

    /// Build vector of mass axes for coordinates in phase space.
    /// Currently only supports two-particle masses.
    /// if argument is left empty, a default set of axes is constructed
    /// \return MassAxes for requested particle combinations
    /// \param pcs vector of vectors of particle indices
    /// \todo Find way to make const.
    const MassAxes massAxes(std::vector<std::vector<unsigned> > pcs = {}) const;

    /// @}

    /// create an empty data set
    /// \param n Number of empty data points to place inside data set
    DataSet createDataSet(size_t n = 0);

    /// Set VariableStatus'es of all Parameter's to unchanged, or leave as fixed
    void setParameterFlagsToUnchanged();

    /// grant friend status to DataAccessor to register itself with this
    friend class DataAccessor;

    /// grant friend status to DataAccessor to register itself with this
    friend class RecalculableDataAccessor;

    /// grant friend status to DataAccessor to register itself with this
    friend class StaticDataAccessor;

    /// grant friend status to DecayingParticle to call addParticleCombination
    friend class DecayingParticle;

protected:

    /// add ParticleCombination to to FourMomenta_
    /// (along with it's daughters through recursive calling) if it is NOT for a FSP.
    virtual void addParticleCombination(const ParticleCombination& pc);

private:

    /// stores whether model structure can be modified
    /// (whether DataAccessors can be added or not)
    bool Locked_;

    /// Lab coordinate system to use in calculating helicity angles
    CoordinateSystem<double, 3> CoordinateSystem_;

    /// ParticleCombination cache
    mutable ParticleCombinationCache ParticleCombinationCache_;

    /// SpinAmplitude cache
    std::unique_ptr<SpinAmplitudeCache> SpinAmplitudeCache_;

    /// Set of all DataAccessor's registered to this model
    DataAccessorSet DataAccessors_;

    /// Vector of raw pointers to all StaticDataAccessor's registered
    /// to this model, in order in which they need be calculated
    StaticDataAccessorVector StaticDataAccessors_;

    /// set of pointers to RecalculableDataAccessors
    RecalculableDataAccessorSet RecalculableDataAccessors_;

    /// Components of full intensity
    std::vector<ModelComponent> Components_;
    
    /// Initial states
    std::vector<std::shared_ptr<DecayingParticle> > InitialStates_;
    
    /// vector of final state particles
    FinalStateParticleVector FinalStateParticles_;

    /// four momenta manager
    std::shared_ptr<FourMomenta> FourMomenta_;

    /// helicity angles manager
    HelicityAngles HelicityAngles_;

};

/// Fix all FreeAmplitude's in a model that parameterize the only
/// possible decay of their parent state.
void fix_solitary_free_amplitudes(Model& m);

/// \return vector of shared_ptr to DecayingParticles inside Model
/// that decay to its full final state, sorted such that the first
/// entries have fixed prefactors
std::vector<std::shared_ptr<DecayingParticle> > full_final_state_isp(const Model& M);

/// \return intensity for all spin projections of an ISP
const double intensity(const ModelComponent& c, const DataPoint& d);

/// \return intensity for a data point evaluated over isp_map
const double intensity(const std::vector<ModelComponent>& C, const DataPoint& d);

/// \return intensity for a data point evaluated over isp_map of model
inline const double intensity(const Model& M, const DataPoint& d)
{ return intensity(M.components(), d); }

/// \return The sum of the logs of squared amplitudes evaluated over the data partition
/// \param M Model to evaluate
/// \param D DataPartition to evalue over
/// \param ped Pedestal to substract from each term in the sum
const double sum_of_log_intensity(const Model& M, DataPartition& D, double ped = 0);

/// \return The sum of the logs of squared amplitudes evaluated over the data partitions
/// \param DP DataPartitionVector of partitions to use
/// \param ped Pedestal to substract from each term in the sum
const double sum_of_log_intensity(const Model& M, DataPartitionVector& DP, double ped = 0);

/// \return all free amplitudes in a model
FreeAmplitudeSet free_amplitudes(const Model& M);

/// \return free amplitude in a model from decay trees evaluating to true
template <typename Last, typename ... UnaryPredicates>
FreeAmplitudeSet free_amplitudes(const Model& M, Last p, UnaryPredicates ... P)
{ return filter(free_amplitudes(M), p, P...); }

/// \return lone free amplitude in a model passing predicates
/// throws if unique amplitude is not found
template <typename Last, typename ... UnaryPredicates>
FreeAmplitudeSet::value_type free_amplitude(const Model& M, Last p, UnaryPredicates ... P)
{ return lone_elt(filter(free_amplitudes(M), p, P...)); }

/// \return set of all particles in model
ParticleSet particles(const Model& M);

/// \return Set of particles in model for which all predicates evaluate true
template <typename Last, typename ... UnaryPredicates>
ParticleSet particles(const Model& M, Last p, UnaryPredicates ... P)
{ return filter(particles(M), p, P...); }

/// \return lone particle in model for which all predicates evaluate true
/// throws if unique particle is not found
template <typename Last, typename ... UnaryPredicates>
ParticleSet::value_type particle(const Model& M, Last p, UnaryPredicates ... P)
{ return lone_elt(filter(particles(M), p, P...)); }

}

#endif
