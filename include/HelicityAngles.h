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

#ifndef yap_HelicityAngles_h
#define yap_HelicityAngles_h

#include "CachedDataValue.h"
#include "StaticDataAccessor.h"

#include <TLorentzVector.h>

namespace yap {

/// \class HelicityAngles
/// \brief Calculates, stores and gives access to helicity angles
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude

class HelicityAngles : public StaticDataAccessor
{
public:

    /// Constructor
    HelicityAngles();

    /// Calculate helicity angles for all possible symmetrization indices
    void calculate(DataPoint& d);

    /// add symmetrizationIndex to SymmetrizationIndices_
    virtual void addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c);

    /// Access helicity angles (const)
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return helicity angles of
    double phi(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
    { return HelicityAngles_.value(0, d, symmetrizationIndex(pc)); }

    /// Access helicity angles (const)
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return helicity angles of
    double theta(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
    { return HelicityAngles_.value(1, d, symmetrizationIndex(pc)); }

protected:

    /// Caclulate Lorentz-transformation for helicity frame
    TLorentzRotation hfTransform(const TLorentzVector& daughterLv);

    /// Transform daughters to helicity frame and calculate helicity angles
    /// Calls this funciton recursively
    void transformDaughters(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, std::vector<TLorentzVector> finalStatesHf);

    /// Helicity angles phi and theta
    CachedDataValue HelicityAngles_;

};

}

#endif
