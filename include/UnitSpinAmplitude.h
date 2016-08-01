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

#ifndef yap_UnitSpinAmplitude_h
#define yap_UnitSpinAmplitude_h

#include "ClebschGordan.h"
#include "Constants.h"
#include "ParticleCombination.h"
#include "Spin.h"
#include "SpinAmplitude.h"

namespace yap {

/// \class UnitSpinAmplitude
/// \brief Implements a spin amplitude that always returns Complex_1
/// \author Daniel Greenwald
/// \ingroup SpinAmplitude
/// \todo Remove actual caching
class UnitSpinAmplitude : public SpinAmplitude
{
public:

    /// \return Complex_1 regardless of input
    const std::complex<double> calc(int two_M, const SpinProjectionVector& two_m,
                                    const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const override
    { return Complex_1; }

    /// \return unit amplitude
    /// \param d dummy DataPoint
    /// \param pc dummy ParticleCombination
    /// \param two_M 2 * spin projection of parent
    /// \param two_m SpinProjectionVector of daughters
    virtual const std::complex<double> amplitude(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                         int two_M, const SpinProjectionVector& two_m) const override
    { return Complex_1; }

    /// \return "unit-valued"
    virtual std::string formalism() const override
    { return "unit-valued"; }

    /// grant SpinAmplitudeCache friend status to call constructor
    friend class SpinAmplitudeCache;

protected:

    /// Constructor
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector of daughters
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin angular momentum
    UnitSpinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s)
        : SpinAmplitude(two_J, two_j, l, two_s, equal_always)
    {
        // check if all spins are 0 for non-res decays
        if (two_j.size() > 2 and (two_J != 0 or l != 0 or two_s != 0 or std::any_of(two_j.begin(), two_j.end(), [](unsigned j){return j != 0;})))
            throw exceptions::Exception("More than 2 daughters, but not all spins are 0", "UnitSpinAmplitude::UnitSpinAmplitude");
        
        // for 2 daughters, loop over all possibilities
        for (auto two_M : projections(two_J)) {
            for (const auto& two_m : projections(two_j)) {
                if (two_m.size() > 2) {
                    addAmplitude(two_M, two_m, true);
                } else {
                    try {
                        if (ClebschGordan::nonzeroCoupling(two_j[0], two_m[0], two_j[1], two_m[1], l, two_s, two_J))
                            addAmplitude(two_M, two_m, true);
                    } catch (const exceptions::InconsistentSpinProjection&) { /* ignore */ }
                }
            }
        }
    }

};

}

#endif
