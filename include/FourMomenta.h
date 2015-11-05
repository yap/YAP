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
    void calculate(DataPoint& d);

    /// Access 4-momentum (const)
    /// \param d DataPoint to get data from
    /// \param i Symmetrization index to access
    const TLorentzVector& p(const DataPoint& d, unsigned i) const
    { return d.FourMomenta_.at(i); }

    /// Access 4-momenutm (const)
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return 4-momentum of
    const TLorentzVector& p(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc)
    { return p(d, symmetrizationIndex(pc)); }

    /// Access invariant mass squared
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return squared mass of
    double m2(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc)
    { return pow(m(d, pc), 2); }

    /// Access invariant mass
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return mass of
    double m(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc)
    {
        if (pc->isFinalStateParticle())
            return FinalStateParticleM_[pc->indices()[0]]->value();
        return M_.value(d, symmetrizationIndex(pc));
    }

    /// Access initial-state 4-momentum (const)
    /// \param d DataPoint to get data from
    const TLorentzVector& initialStateMomentum(const DataPoint& d)
    { return p(d, InitialStateIndex_); }

protected:

    /// Symmetrization index of initial state
    int InitialStateIndex_;

    /// mass [GeV]
    RealCachedDataValue M_;

    /// masses of the final state particles
    std::vector<std::shared_ptr<RealParameter> > FinalStateParticleM_;

};

}
#endif

