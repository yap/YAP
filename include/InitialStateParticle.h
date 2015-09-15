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

    /// Constructor
    InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

    ~InitialStateParticle()
    { DataAccessors_.clear(); }

    //virtual Amp amplitude(DataPoint& d) override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Check consistency of data point with this initial state and its final state particles
    /// \param d DataPoint to check
    /// \todo Flesh out
    virtual bool consisent(const DataPoint& d) const
    { return true; }

    /// you MUST call this function after you have added all decay channels and before adding DataPoints
    bool prepare();

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

    /// @}

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

    /// Set parents of symmetrization indices (recursively)
    virtual void setSymmetrizationIndexParents() override;

    /// add DataAccessor to set
    void addDataAccessor(DataAccessor* d)
    { DataAccessors_.insert(d); }

    /// add DataAccessor to set
    void removeDataAccessor(DataAccessor* d)
    { DataAccessors_.erase(d); }

    void setDataAcessorIndices();

    friend class DataAccessor;

    /// List of all DataAccessor objects in the InitialsStateParticle and below
    std::set<DataAccessor*> DataAccessors_;

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
