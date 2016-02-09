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
#include "FourVector.h"
#include "ParticleCombination.h"
#include "StaticDataAccessor.h"

#include <ostream>

namespace yap {

class DataPoint;
class InitialStateParticle;

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
    FourMomenta(InitialStateParticle* isp);

    /// check consistency
    bool consistent() const;

    /// Fill 4-momenta
    /// \param d DataPoint to fill
    /// \param dataPartitionIndex for status tracking
    virtual void calculate(DataPoint& d, unsigned dataPartitionIndex = 0) override;

    /// Access 4-momenutm (const)
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return 4-momentum of
    const FourVector<double>& p(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// Access invariant mass squared
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return squared mass of
    double m2(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
    { return pow(m(d, pc), 2); }

    /// Access invariant mass
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return mass of
    double m(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// Access initial-state 4-momentum (const)
    /// \param d DataPoint to get data from
    const FourVector<double>& initialStateMomentum(const DataPoint& d)
    { return p(d, InitialStatePC_); }

    /// \return masses
    std::shared_ptr<RealCachedDataValue> masses()
    { return M_; }

    /// \return masses (const)
    std::shared_ptr<RealCachedDataValue> masses() const
    { return M_; }

    /// print all masses
    std::ostream& printMasses(const DataPoint& d, std::ostream& os = std::cout) const;

    virtual std::string data_accessor_type() const override
    {return "FourMomenta"; }

    /// grant friend status to InitialStateParticle to call setFourMomenta
    friend class InitialStateParticle;

protected:

    /// looks for ISP when adding ParticleCombination's
    void addParticleCombination(std::shared_ptr<ParticleCombination> pc) override;

    /// override to do nothing, since FourMomenta doesn't rely on parents being set.
    void pruneSymmetrizationIndices() override
    {}

private:

    /// Symmetrization index of initial state
    std::shared_ptr<ParticleCombination> InitialStatePC_;

    /// invariant mass of particle combinations [GeV]
    std::shared_ptr<RealCachedDataValue> M_;

};

}
#endif

