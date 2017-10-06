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

#ifndef yap_Particle_h
#define yap_Particle_h

#include "fwd/Particle.h"

#include "fwd/DataPoint.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/Model.h"
#include "fwd/Parameter.h"
#include "fwd/Particle.h"
#include "fwd/ParticleTable.h"
#include "fwd/Spin.h"

#include "ParticleCombination.h"
#include "QuantumNumbers.h"

#include <iostream>
#include <memory>
#include <string>

namespace yap {

/// \class Particle
/// \brief Abstract Particle base class.
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Particle Particle-related classes
class Particle : public std::enable_shared_from_this<Particle>
{
protected:

    /// Constructor
    /// \param name Name of particle
    /// \param q Quantum numbers of particle
    Particle(const std::string& name, const QuantumNumbers& q)
        : std::enable_shared_from_this<Particle>(),
          QuantumNumbers_(q), Name_(name) {}

    /// Constructor
    /// \param pde ParticleTableEntry to take name and quantum numbers from
    Particle(const ParticleTableEntry& pde);
    
public:

    /// Check consitency of object
    virtual bool consistent() const;

    /// \name Getters
    /// @{

    /// const access QuantumNumbers object
    const QuantumNumbers& quantumNumbers() const
    { return QuantumNumbers_; }

    /// Get name
    const std::string& name() const
    { return Name_; }

    /// get raw pointer to Model (const)
    virtual const Model* model() const = 0;

    /// \return ParticleCombinations_
    const ParticleCombinationSet& particleCombinations() const
    { return ParticleCombinations_; }

    /// @}

    /// grant friend status to DecayChannel to call addParticleCombination
    friend class DecayChannel;

protected:

    /// add ParticleCombination to ParticleCombinations_
    virtual void addParticleCombination(const ParticleCombination& pc) = 0;

    /// prune ParticleCombinations_ to only contain ParticleCombination's tracing back up the ISP
    virtual void prune()
    { prune_particle_combinations(ParticleCombinations_); }

    /// register any necessary DataAccessor's with model
    virtual void registerWithModel() = 0;

private:

    /// Quantum numbers of particle
    QuantumNumbers QuantumNumbers_;

    /// Name of particle
    std::string Name_;

    /// set of ParticleCombinations that can comprise this particle
    ParticleCombinationSet ParticleCombinations_;

};

/// \return charge of ParticleVector
inline const int charge(const ParticleVector& P)
{ return std::accumulate(P.begin(), P.end(), 0, [](int q, const auto& p){return q + p->quantumNumbers().Q();}); }

/// \return parity of ParticleVector
inline const int parity(const ParticleVector& P)
{ return std::accumulate(P.begin(), P.end(), 1, [](int q, const auto& p){return q * p->quantumNumbers().P();}); }

/// \return SpinVector from ParticleVector
const SpinVector spins(const ParticleVector& v);

/// \return whether a particle decays to its model's full final state
const bool decays_to_full_final_state(const Particle& p);

/// convert to string
inline std::string to_string(const Particle& p)
{ return p.name() + "(" + to_string(p.quantumNumbers()) + ")"; }

/// convert to string
inline std::string to_string(const ParticleVector& p)
{ return std::accumulate(p.begin(), p.end(), std::string(""), [](std::string& s, const auto& p){return s += ", " + p->name();}).erase(0, 2); }

/// convert to string
std::string to_string(const ParticleVector& p, const SpinProjectionVector& two_m);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const Particle& p)
{ os << to_string(p); return os; }

}

#endif
