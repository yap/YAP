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

#include "SpinAmplitude.h"

#include <complex>
#include <map>
#include <memory>

namespace yap {

class InitialStateParticle;
class ParticleCombination;

/// \class HelicitySpinAmplitude
/// \brief Class implementing a canonical spin amplitude, i.e. with defined relative angular momentum.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude

class HelicitySpinAmplitude : public SpinAmplitude
{
public:

    /// \name Constructors
    /// @{

    /// Constructor
    HelicitySpinAmplitude(const QuantumNumbers& initial,
                          const QuantumNumbers& final1, const QuantumNumbers& final2, unsigned char twoL);

    /// @}

    /// Calculate complex amplitude
    virtual std::complex<double> amplitude(DataPartition& d, const std::shared_ptr<const ParticleCombination>& pc) const override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// also calculates ClebschGordan coefficient
    virtual void addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c);

    /// also clears ClebschGordanCoefficients_
    virtual void clearSymmetrizationIndices();

    /// cast into string
    operator std::string() const override;

    /// Calculate Clebsch-Gordan coefficients for all particleCombinations
    double calculateClebschGordanCoefficient(std::shared_ptr<const ParticleCombination> c) const;

    // virtual std::vector<std::shared_ptr<ComplexParameter> > ParametersItDependsOn() override;

    virtual CachedDataValueSet CachedDataValuesItDependsOn() override
    { return {SpinAmplitude_}; }


private:

    /// Check if SpinAmplitudes are equal
    bool equals(const SpinAmplitude& rhs) const override;

    /// Clebsch-Gordan coefficient for 2*λ_1, 2*λ_2
    /// \todo make this a Parameter???
    std::map<std::shared_ptr<const ParticleCombination>, double> ClebschGordanCoefficients_;

    std::shared_ptr<ComplexCachedDataValue> SpinAmplitude_;

};

}

#endif
