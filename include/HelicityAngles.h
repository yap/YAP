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
#include "FourVector.h"
#include "LorentzTransformation.h"
#include "StaticDataAccessor.h"

namespace yap {

/// \class HelicityAngles
/// \brief Calculates, stores and gives access to helicity angles
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude
///
/// Angles are defined as follows:
///   -# Given decaying particle with momentum \f$ \vec{p} \f$ in its parent frame and
///      the \f$ z \f$ axis of its parent frame \f$ \hat{z}_0 \f$,
///      define the reference frame:
///      - \f$ \hat{z} \equiv \hat{p} \f$
///      - \f$ \hat{y} \equiv \hat{z}_0 \cross \hat{z} \f$
///      - \f$ \hat{x} \equiv \hat{y} \cross \hat{z} \f$
///   -# Given daughter particle with momentum \f$ \vec{q} \f$
///      in the decaying particle's rest frame, the angles are
///      - \f$ \cos\theta \equiv \hat{q} \cdot \hat{z} \f$
///      - \f$ \cos\phi   \equiv \hat{q} \cdot \hat{x} / \sin\theta \f$
/// If the decaying particle has 0 momentum or \f$ \hat{p} = \hat{z}_0 \f$
///   - \f$ \hat{z} \equiv \hat{z}_0 \f$
///   - \f$ \hat{y} \equiv \hat{y}_0 \f$
///   - \f$ \hat{x} \equiv \hat{x}_0 \f$
/// with the 0th coordinate system given by the user for the inital state,
/// and defaulting to standard Cartesian system as defined in #Constants.h

class HelicityAngles : public StaticDataAccessor
{
public:

    /// Constructor
    HelicityAngles();

    /// Calculate helicity angles for all possible symmetrization indices
    virtual void calculate(DataPoint& d) override;

    // /// add symmetrizationIndex to SymmetrizationIndices_
    // virtual void addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c);

    /// get azimuthal angle
    double phi(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
    { return Phi_->value(d, symmetrizationIndex(pc)); }

    /// access azimuthal angle
    std::shared_ptr<RealCachedDataValue>& phi()
    { return Phi_; }

    /// access azimuthal angle (const)
    const std::shared_ptr<RealCachedDataValue>& phi() const
    { return Phi_; }

    /// get polar angle
    double theta(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
    { return Theta_->value(d, symmetrizationIndex(pc)); }

    /// access polar angle
    std::shared_ptr<RealCachedDataValue>& theta()
    { return Theta_; }

    /// access polar angle (const)
    const std::shared_ptr<RealCachedDataValue>& theta() const
    { return Theta_; }

protected:

    /// recursive helicity-angle calculator that travels down decay trees for all channels
    void calculateAngles(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts);

    /// Azimuthal angle
    std::shared_ptr<RealCachedDataValue> Phi_;

    /// Polar angle
    std::shared_ptr<RealCachedDataValue> Theta_;

};

}

#endif
