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

namespace yap {

class DataSet;

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

    /// Clone
    virtual std::shared_ptr<Particle> clone() const override
    { return std::make_shared<InitialStateParticle>(*this); }

    /// @}

    /// \todo remove!
    double logLikelihood();

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Check consistency of data point with this initial state and its final state particles
    /// \param d DataPoint to check
    /// \todo Flesh out
    virtual bool consisent(const DataPoint& d) const
    { return true; }

    /// you MUST call this function after you have added all decay channels and before adding DataPoints
    bool prepare();

    /// set free amplitudes to DecayChannels
    bool setFreeAmplitudes(const std::vector<Amp>& amps);

    /// \name Getters
    /// @{

    FourMomenta& fourMomenta()
    { return FourMomenta_; }

    const FourMomenta& fourMomenta() const
    { return FourMomenta_; }

    HelicityAngles& helicityAngles()
    { return HelicityAngles_; }

    const HelicityAngles& helicityAngles() const
    { return HelicityAngles_; }

    bool prepared() const
    {return Prepared_; }

    /// \return free amplitudes of DecayChannels_
    std::vector<Amp> freeAmplitudes() const;

    DataSet& dataSet()
    { return DataSet_; }

    /// @}

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

    void printDataAccessors();

private:

    friend void DecayingParticle::optimizeSpinAmplitudeSharing();

    /// Set parents of symmetrization indices (recursively)
    virtual void setSymmetrizationIndexParents() override;

    /// add DataAccessor to set
    void addDataAccessor(DataAccessor* d)
    { DataAccessors_.insert(d); }

    /// remove DataAccessor from set
    void removeDataAccessor(DataAccessor* d);

    void setDataAcessorIndices();

    friend class DataAccessor;

    bool Prepared_;

    /// List of all DataAccessor objects in the InitialsStateParticle and below
    std::set<DataAccessor*> DataAccessors_;

    /// List of all DecayChannel objects in the InitialsStateParticle and below
    std::vector<DecayChannel*> DecayChannels_;

    /// vector of final state particles
    std::vector<FinalStateParticle*> FinalStateParticles;

    /// four momenta manager
    FourMomenta FourMomenta_;

    /// helicity angles manager
    HelicityAngles HelicityAngles_;

    /// Data set
    DataSet DataSet_;

};

}

#endif
