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

#ifndef yap_FourMomenta_
#define yap_FourMomenta_

#include "CachedDataValue.h"
#include "ParticleCombination.h"
#include "StaticDataAccessor.h"

#include <TLorentzVector.h>

namespace yap {

class DataPoint;

/// \class FourMomenta
/// \brief Stores and gives access to four-momenta and invariant masses
/// \author Johannes Rauch, Daniel Greenwald
///
/// The final-state particles must have indices corresponding to number;
/// initial state must also exist in entries

class FourMomenta : public StaticDataAccessor
{
public:

    /// Constructor
    FourMomenta();

    /// Find ISP in set and store index location
    /// Fill FinalStateParticleM_ and FinalStateParticleM2_
    void prepare();

    /// check consistency
    bool consistent() const;

    /// Fill 4-momenta
    virtual void calculate(DataPoint& d) override;

    /// Access 4-momenutm (const)
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return 4-momentum of
    const TLorentzVector& p(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc)
    {
        if (pc->isFinalStateParticle())
            return d.FSPFourMomenta_[pc->indices()[0]];
        return d.FourMomenta_[symmetrizationIndex(pc)];
    }

    /// Access invariant mass squared
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return squared mass of
    double m2(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
    { return pow(m(d, pc), 2); }

    /// Access invariant mass
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return mass of
    double m(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
    {
        if (pc->isFinalStateParticle())
            return FinalStateParticleM_[pc->indices()[0]]->value();
        return M_->value(d, symmetrizationIndex(pc));
    }

    /// Access initial-state 4-momentum (const)
    /// \param d DataPoint to get data from
    const TLorentzVector& initialStateMomentum(const DataPoint& d)
    { return p(d, InitialStatePC_); }

    /// \name Dalitz coordinate stuff
    /// @{

    /// get ParticleCombinationVector for requested combinations
    ParticleCombinationVector getDalitzAxes(std::vector<std::vector<ParticleIndex> > pcs) const;

    /// get pair masses
    ParticleCombinationMap<double> pairMasses(const DataPoint& d) const;

    /// get pair masses squared
    ParticleCombinationMap<double> pairMassSquares(const DataPoint& d) const;

    /// set masses of particle combinations in axes to those in masses
    /// \param d DataPoint to set into
    /// \param axes vector of ParticleCombination's to set masses of
    /// \param masses vector of masses to be set to
    /// \return success of action
    bool setMasses(DataPoint& d, const ParticleCombinationVector& axes, const std::vector<double>& masses);

    /// set masses of particle combinations in axes to those in masses
    /// \param d DataPoint to set into
    /// \param axes vector of ParticleCombination's to set masses of
    /// \param masses vector of masses to be set to
    /// \return success of action
    bool setSquaredMasses(DataPoint& d, const ParticleCombinationVector& axes, const std::vector<double>& squaredMasses);

    /// set masses
    /// also calculate remaining masses and final-state 4-momenta
    /// \return if masses form a complete set and are in phase space
    bool setMasses(DataPoint& d, ParticleCombinationMap<double> m);

    /// set masses squared
    /// also calculate remaining masses and final-state 4-momenta
    /// \return if masses form a complete set and are in phase space
    bool setMassSquares(DataPoint& d, ParticleCombinationMap<double> m2);

    /// @}

    std::shared_ptr<RealCachedDataValue> masses()
    { return M_; }

    /// print all masses
    void printMasses(const DataPoint& d) const;

protected:

    /// set all masses to -1 (except FinalStateParticleM_)
    void resetMasses(DataPoint& d);

    /// calculate all masses from a complete set of masses
    /// \return success of action
    bool calculateMissingMasses(DataPoint& d);

    /// calculate four-momenta from squared invariant masses
    /// with the following convention for three-momenta:\n
    /// p1 defines +z direction
    /// p1 x p2 defines +y direction
    std::vector<TLorentzVector> calculateFourMomenta(const DataPoint& d) const;

    /// \return set of all pair particle combinations, without duplicates
    ParticleCombinationVector pairParticleCombinations() const;

    /// Symmetrization index of initial state
    std::shared_ptr<const ParticleCombination> InitialStatePC_;

    /// \todo Perhaps find better way than storing FinalStatePC, RecoilPC, and PairPC.

    /// Final-state particle PCs
    ParticleCombinationVector FinalStatePC_;

    /// Symmetrization indices for recoil against each final-state particle,
    /// index in vector is FSP index
    ParticleCombinationVector RecoilPC_;

    /// Symmetrization indices for pairs of particles
    std::vector<ParticleCombinationVector> PairPC_;

    /// invariant mass of particle combinations [GeV]
    std::shared_ptr<RealCachedDataValue> M_;

    /// masses of the final state particles
    std::vector<std::shared_ptr<RealParameter> > FinalStateParticleM_;

};

}
#endif

