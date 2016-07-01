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

#include "fwd/DataAccessor.h"
#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/DataSet.h"
#include "fwd/DecayingParticle.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/FourMomenta.h"
#include "fwd/FourVector.h"
#include "fwd/HelicityAngles.h"
#include "fwd/MassAxes.h"
#include "fwd/MeasuredBreakupMomenta.h"
#include "fwd/Parameter.h"
#include "fwd/RecalculableDataAccessor.h"
#include "fwd/SpinAmplitudeCache.h"
#include "fwd/StaticDataAccessor.h"
#include "fwd/StatusManager.h"

#include "CoordinateSystem.h"
#include "ParticleCombinationCache.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

/// map initial state particle to free (real) amplitude (for incoherent summing over initial state particles)
using initialStateParticleMap = std::map<std::shared_ptr<DecayingParticle>, std::shared_ptr<RealParameter> >;

/// \class Model
/// \brief Class implementing a PWA model
/// \author Johannes Rauch, Daniel Greenwald
class Model
{
public:

    /// Constructor
    /// \param SAC unique_ptr to SpinAmplitudeCache
    Model(std::unique_ptr<SpinAmplitudeCache> SAC);

    /// Calculate model for each data point in the data partition
    /// \param D DataPartition to calculate over
    /// \todo This need not be a member function!
    void calculate(DataPartition& D) const;

    /// Check consistency of object
    virtual bool consistent() const;

    bool locked() const
    { return Locked_; }

    /// removes expired DataAccessor's, prune's remaining, and assigns them indices
    void prepareDataAccessors();

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

    /// \return MeasuredBreakupMomenta accessor
    std::shared_ptr<MeasuredBreakupMomenta> measuredBreakupMomenta()
    { return MeasuredBreakupMomenta_; }

    /// \return MeasuredBreakupMomenta accessor (const)
    const std::shared_ptr<MeasuredBreakupMomenta> measuredBreakupMomenta() const
    { return MeasuredBreakupMomenta_; }

    /// \return HelicityAngles accessor
    /// is nullptr if not using helicity formalism
    std::shared_ptr<HelicityAngles> helicityAngles()
    { return HelicityAngles_; }

    /// \return HelicityAngles accessor (const)
    /// is nullptr if not using helicity formalism
    const std::shared_ptr<HelicityAngles> helicityAngles() const
    { return HelicityAngles_; }

    /// \return ParticleCombinationCache
    ParticleCombinationCache& particleCombinationCache()
    { return ParticleCombinationCache_; }

    /// \return ParticleCombinationCache (const)
    const ParticleCombinationCache& particleCombinationCache() const
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

    const initialStateParticleMap& initialStateParticles() const
    { return InitialStateParticles_; }

    /// \return THE initial state particle
    /// Per default, this will be the first initial state particle added, that decays to the full final state
    const std::shared_ptr<DecayingParticle> initialStateParticle() const
    { return TheInitialStateParticle_; }

    /// \return set of DataAccessors
    const DataAccessorSet& dataAccessors() const
    { return DataAccessors_; }

    /// \return set of DataAccessors
    const StaticDataAccessorVector& staticDataAccessors() const
    { return StaticDataAccessors_; }

    /// \return (min, max) array[2] of mass range for particle combination
    /// \param pc shared pointer to ParticleCombination to get mass range of
    /// \param p Initial state particle to get the mass range for
    std::array<double, 2> massRange(const std::shared_ptr<ParticleCombination>& pc, std::shared_ptr<DecayingParticle> initialStateParticle) const;

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

    /// add an initial state particle
    /// The first particle added that decays to the full final state will become THE initial state particle,
    /// which can be retrieved by calling #initialStateParticle,
    /// and its amplitude will be fixed
    const initialStateParticleMap::value_type& addInitialStateParticle(std::shared_ptr<DecayingParticle> bg);

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
    const MassAxes massAxes(std::vector<std::vector<unsigned> > pcs = {});

    /// Calculate four-momenta for final-state particles for phase-space coordinate
    /// \param axes phase-space axes
    /// \param squared_masses phase-space coordinate
    /// \param p Initial state particle to calculate the four momenta with (its mass will be used)
    std::vector<FourVector<double> > calculateFourMomenta(const MassAxes& axes, const std::vector<double>& squared_masses, std::shared_ptr<DecayingParticle> initialStateParticle) const;

    /// @}

    /// create an empty data set
    /// \param n Number of empty data points to place inside data set
    DataSet createDataSet(size_t n = 0);

    /// Set VariableStatus'es of all Parameter's to unchanged, or leave as fixed
    void setParameterFlagsToUnchanged();

    /// Print the list of DataAccessor's
    void printDataAccessors(bool printParticleCombinations = true) const;

    /// Print all VariableStatus'es and CalculationStatus'es
    void printFlags(const StatusManager& sm) const;

    /// grant friend status to DataAccessor to register itself with this
    friend class DataAccessor;

    /// grant friend status to DecayingParticle to call addParticleCombination
    friend class DecayingParticle;

protected:

    /// add ParticleCombination to to FourMomenta_, HelicityAngles_, and
    /// MeasuredBreakupMomenta_ (along with it's daughters through
    /// recursive calling) if it is NOT for a FSP.
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> pc);

    /// register a DataAccessor with this Model
    virtual void addDataAccessor(DataAccessorSet::value_type da);

    /* /// remove a DataAccessor from this Model */
    /* virtual void removeDataAccessor(DataAccessorSet::value_type da); */

private:

    /// lock model
    void lock()
    { Locked_ = true; }

    /// stores whether model structure can be modified
    /// (whether DataAccessors can be added or not)
    bool Locked_;

    /// Lab coordinate system to use in calculating helicity angles
    CoordinateSystem<double, 3> CoordinateSystem_;

    /// ParticleCombination cache
    ParticleCombinationCache ParticleCombinationCache_;

    /// SpinAmplitude cache
    std::unique_ptr<SpinAmplitudeCache> SpinAmplitudeCache_;

    /// Set of all DataAccessor's registered to this model
    DataAccessorSet DataAccessors_;

    /// Vector of raw pointers to all StaticDataAccessor's registered
    /// to this model, in order in which they need be calculated
    StaticDataAccessorVector StaticDataAccessors_;

    /// set of pointers to RecalculableDataAccessors
    RecalculableDataAccessorSet RecalculableDataAccessors_;

    /// pointers to initial particles
    /// they will be summed in incoherently
    initialStateParticleMap InitialStateParticles_;

    /// The main initial state particle.
    /// Per default, this will be the first initialStateParticle added, that decays to the full final state
    std::shared_ptr<DecayingParticle> TheInitialStateParticle_;

    /// vector of final state particles
    FinalStateParticleVector FinalStateParticles_;

    /// four momenta manager
    std::shared_ptr<FourMomenta> FourMomenta_;

    /// Breakup momenta manager
    std::shared_ptr<MeasuredBreakupMomenta> MeasuredBreakupMomenta_;

    /// helicity angles manager
    std::shared_ptr<HelicityAngles> HelicityAngles_;

};

/// \return The sum of the logs of squared amplitudes evaluated over the data partition
/// \param M Model to evaluate
/// \param D DataPartition to evalue over
const double sumOfLogsOfSquaredAmplitudes(const Model& M, DataPartition& D);

/// \return The sum of the logs of squared amplitudes evaluated over the data partitions
/// \param DP DataPartitionVector of partitions to use
const double sumOfLogsOfSquaredAmplitudes(const Model& M, DataPartitionVector& DP);

}

#endif
