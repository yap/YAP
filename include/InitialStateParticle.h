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

#ifndef yap_InitialStateParticle_h
#define yap_InitialStateParticle_h

#include "DataPartition.h"
#include "DataSet.h"
#include "DecayingParticle.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "MeasuredBreakupMomenta.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

/// \class InitialStateParticle
/// \brief Class implementing an initial state particle.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class InitialStateParticle : public DecayingParticle
{
public:

    /// \name Constructor, destructor, clone
    /// @{

    /// Constructor
    InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

    /// Destructor
    ~InitialStateParticle();

    /// @}

    /// \name InitialStateParticle friends
    /// @{

    friend class DataAccessor;
    friend class AmplitudeComponent;
    friend void DecayingParticle::optimizeSpinAmplitudeSharing();

    /// @}

    /// \todo remove/rename/rework!
    double logLikelihood(DataPartition& d);

    /// Check consistency of object
    virtual bool consistent() const override;

    /// you MUST call this function after you have added all decay channels and before adding DataPoints
    bool prepare();

    /// \name Getters
    /// @{

    FourMomenta& fourMomenta()
    { return FourMomenta_; }

    const FourMomenta& fourMomenta() const
    { return FourMomenta_; }

    MeasuredBreakupMomenta& measuredBreakupMomenta()
    { return MeasuredBreakupMomenta_; }

    const MeasuredBreakupMomenta& measuredBreakupMomenta() const
    { return MeasuredBreakupMomenta_; }

    HelicityAngles& helicityAngles()
    { return HelicityAngles_; }

    const HelicityAngles& helicityAngles() const
    { return HelicityAngles_; }

    bool prepared() const
    {return Prepared_; }

    /// \return free amplitudes of DecayChannels_
    std::vector<std::shared_ptr<ComplexParameter> > freeAmplitudes() const;

    DataSet& dataSet()
    { return DataSet_; }

    /// @}

    /// \name DataPoints
    /// @{

    /// Add data point via four-momenta
    /// This method is faster since it avoids unneccessary copying of objects
    /// and resizing of the DataPoint's storage
    bool addDataPoint(const std::vector<TLorentzVector>& fourMomenta);

    /// Add data point via move
    /// \param d DataPoint to move into DataSet
    /// \return Success of action
    bool addDataPoint(DataPoint&& d);

    /// Add data point via copy
    /// \param d DataPoint to copy into DataSet
    /// \return Success of action
    bool addDataPoint(const DataPoint& d);

    /// @}

    void printDataAccessors(bool printParticleCombinations = true);

private:

    /// Set parents of symmetrization indices (recursively)
    virtual void setSymmetrizationIndexParents() override;

    /// set all CachedDataValue flags to kNeedsCheck
    /// call before looping over a DataPartition
    void setParameterCachedDataValueFlagsToNeedsCheck(unsigned dataPartitionIndex);

    /// set all parameter flags to kUnchanged (or leave at kFixed)
    /// call after looping over a DataPartition
    void setParameterFlagsToUnchanged(unsigned dataPartitionIndex);

    /// add DataAccessor to set
    void addDataAccessor(DataAccessor* d);

    /// remove DataAccessor from set
    void removeDataAccessor(DataAccessor* d);

    /// add AmplitudeComponent to set
    void addAmplitudeComponent(AmplitudeComponent* d)
    { AmplitudeComponents_.insert(d); }

    /// remove AmplitudeComponent from set
    void removeAmplitudeComponent(AmplitudeComponent* d);

    void setDataAcessorIndices();

    bool Prepared_;

    /// List of all DataAccessor objects in the InitialsStateParticle and below
    std::set<DataAccessor*> DataAccessors_;

    /// List of all AmplitudeComponent objects in the InitialsStateParticle and below
    std::set<AmplitudeComponent*> AmplitudeComponents_;

    /// List of all DecayChannel objects in the InitialsStateParticle and below
    std::vector<DecayChannel*> DecayChannels_;

    /// vector of final state particles
    std::vector<FinalStateParticle*> FinalStateParticles;

    /// four momenta manager
    FourMomenta FourMomenta_;

    /// Breakup momenta manager
    MeasuredBreakupMomenta MeasuredBreakupMomenta_;

    /// helicity angles manager
    HelicityAngles HelicityAngles_;

    /// Data set
    DataSet DataSet_;

};

}

#endif
