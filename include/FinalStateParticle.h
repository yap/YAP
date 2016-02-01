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

#ifndef yap_FinalStateParticle_h
#define yap_FinalStateParticle_h

#include "CalculationStatus.h"
#include "Constants.h"
#include "DataAccessor.h"
#include "DataPoint.h"
#include "InitialStateParticle.h"
#include "Particle.h"
#include "ParticleIndex.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

class ParticleCombination;

/// \class FinalStateParticle
/// \brief Class representing a final-state particle
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class FinalStateParticle : public Particle
{
public:

    /// \name Constructor
    /// @{

    /// Constructor
    /// \param q Quantum numbers of particle
    /// \param m Mass of particle
    /// \param name Name of particle
    FinalStateParticle(const QuantumNumbers& q, double m, std::string name);

    /// @}

    /// Calculate complex amplitude.
    /// All parameters are ignored. See particle::amplitude() for info.
    /// \return 1 + 0i
    virtual std::complex<double> amplitude(DataPoint&, const std::shared_ptr<ParticleCombination>&, int, unsigned) const override
    { return Complex_1; }

    /// Check consistency
    virtual bool consistent() const override;

    /// \return list of all ParticleCombinations
    ParticleCombinationVector particleCombinations() const override
    { return ParticleCombinations_; }

    /// \return owning InitialStateParticle
    InitialStateParticle* initialStateParticle() override
    { return InitialStateParticle_; }

    /// \name Friends
    /// @{

    /// Grant ISP friendship to set FSP's indices
    friend void InitialStateParticle::setFinalStateParticles(std::initializer_list<std::shared_ptr<FinalStateParticle> >);

    /// @}

protected:

    /// set raw pointer to owning InitialStateParticle
    void setInitialStateParticle(InitialStateParticle* isp)
    { InitialStateParticle_ = isp; }

    /// add ParticleCombination to ParticleCombinations
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> pc) override;

    /// \return empty set
    virtual DataAccessorSet dataAccessors() override
    { return DataAccessorSet(); }

private:

    /// raw pointer to initial state particle decaying to this final state particle
    InitialStateParticle* InitialStateParticle_;

    ParticleCombinationVector ParticleCombinations_;

};

}

#endif
