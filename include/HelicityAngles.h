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

#include "fwd/CachedValue.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"

#include "Constants.h"
#include "CoordinateSystem.h"
#include "Rotation.h"
#include "StaticDataAccessor.h"
#include "ThreeVector.h"

#include <memory>

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
///      - \f$ \hat{y} \equiv \hat{z}_0 \times \hat{z} \f$
///      - \f$ \hat{x} \equiv \hat{y} \times \hat{z} \f$
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
    /// \param m owning Model
    HelicityAngles(Model& m);

    /// Calculate helicity angles for all possible symmetrization indices
    /// \param d DataPoint to calculate into
    /// \param sm StatusManager to update
    virtual void calculate(DataPoint& d, StatusManager& sm) const override;

    /// get azimuthal angle
    double phi(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const;

    /// access azimuthal angle
    std::shared_ptr<RealCachedValue>& phi()
    { return Phi_; }

    /// access azimuthal angle (const)
    const std::shared_ptr<RealCachedValue>& phi() const
    { return Phi_; }

    /// get polar angle
    double theta(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const;

    /// access polar angle
    std::shared_ptr<RealCachedValue>& theta()
    { return Theta_; }

    /// access polar angle (const)
    const std::shared_ptr<RealCachedValue>& theta() const
    { return Theta_; }

    /// grant friend status to Model to call addParticleCombination
    friend class Model;

    /// grant friend status to DataAccessor to call addParticleCombination
    friend class DataAccessor;

protected:

    /// add to model's StaticDataAccessors_
    void virtual addToStaticDataAccessors() override;

    /// recursive helicity-angle calculator that travels down decay trees for all channels
    void calculateAngles(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc,
                         const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts,
                         StatusManager& sm) const;

    /// override to throw on adding non-two-body PC
    void addParticleCombination(const ParticleCombination& pc) override;

private:

    /// Azimuthal angle
    std::shared_ptr<RealCachedValue> Phi_;

    /// Polar angle
    std::shared_ptr<RealCachedValue> Theta_;

};

/// Calculate helicity frame of V transformed from C,
/// with z = unit(V), y = C.z X z, x = y X z
/// \param V ThreeVector defining new Z direction
/// \param C CoordinateSystem aiding in defining new Y direction
template <typename T>
CoordinateSystem<T, 3> helicityFrame(const ThreeVector<T>& V, const CoordinateSystem<T, 3>& C)
{
    constexpr T epsilon = T(5) * std::numeric_limits<T>::epsilon();

    // if V is 0, return C
    if (norm(V) <= epsilon)
        return C;

    CoordinateSystem<double, 3> vC;

    // Set Z to V direction
    vC[2] = unit(V);

    // calc angle between new Z and old Z
    auto theta = angle(vC[2], C[2]);

    // if Z direction is same, return C
    if (std::abs(theta) <= epsilon)
        return C;

    // if Z direction is opposite, return C rotated 180 degrees around Y axis
    if (std::abs(theta - pi<double>()) <= epsilon)
        return rotation<T>(C[1], rad(180.)) * C;

    // Y := (C's Z) cross Z
    vC[1] = unit(cross(C[2], vC[2]));

    // X := Y cross Z (right-handed)
    vC[0] = cross(vC[1], vC[2]);

    return vC;
}


}

#endif
