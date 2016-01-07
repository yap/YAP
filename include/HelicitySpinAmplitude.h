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

#ifndef yap_HelicitySpinAmplitude_h
#define yap_HelicitySpinAmplitude_h

#include "DataPoint.h"
#include "ParticleCombination.h"
#include "QuantumNumbers.h"
#include "SpinAmplitude.h"

#include <complex>
#include <map>
#include <memory>

namespace yap {

class ParticleCombination;

template <class T>
class SpinAmplitudeCache;

/// \class HelicitySpinAmplitude
/// \brief Class implementing a canonical spin amplitude, i.e. with defined relative angular momentum.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude
class HelicitySpinAmplitude : public SpinAmplitude
{
public:

    /// Calculate spin amplitude for all possible symmetrization indices
    virtual void calculate(DataPoint& d) override;

    /// Check consistency of object
    virtual bool consistent() const
    { return true; }

    /// check if Clebsch-Gordan coefficient is nonzero before adding pc,
    /// also add all helicity states of the parent
    virtual ParticleCombinationVector addSymmetrizationIndices(std::shared_ptr<ParticleCombination> pc) override;

    /// also calculates ClebschGordan coefficient
    virtual void addSymmetrizationIndex(std::shared_ptr<ParticleCombination> pc);

    /// also clears ClebschGordanCoefficients_
    virtual void clearSymmetrizationIndices();

    /// cast into string
    operator std::string() const override
    { return SpinAmplitude::operator std::string() + " in helicity formalism"; }

    /// grant SpinAmplitudeCache friend status to call constructor
    friend class SpinAmplitudeCache<HelicitySpinAmplitude>;

protected:

    /// Constructor
    /// \param intial quantum numbers of Initial-state
    /// \param final1 quantum numbers of first daughter
    /// \param final2 quantum numbers of second daughter
    /// \param orbital angular momentum
    /// \param isp raw pointer to owning InitialStateParticle
    HelicitySpinAmplitude(const QuantumNumbers& initial,
                          const QuantumNumbers& final1,
                          const QuantumNumbers& final2,
                          unsigned l,
                          InitialStateParticle* isp);

private:
    /// check equality
    virtual bool equals(const SpinAmplitude& other) const override;

    /// spin for l-s coupling.
    /// \todo Why is j1 + j2?
    unsigned two_s()
    { return finalQuantumNumbers()[0].twoJ() + finalQuantumNumbers()[1].twoJ(); }

    /// Clebsch-Gordan coefficient for 2*λ_1, 2*λ_2
    ParticleCombinationMap<double> ClebschGordanCoefficients_;

    std::shared_ptr<ComplexCachedDataValue> SpinAmplitude_;

};

}

#endif
