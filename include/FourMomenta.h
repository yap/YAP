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

#include "FourVector.h"
#include "StaticDataAccessor.h"

#include <memory>
#include <ostream>
#include <vector>

namespace yap {

class DataPoint;
class FourVectorCachedDataValue;
class Model;
class ParticleCombination;
class RealCachedDataValue;

/// \class FourMomenta
/// \brief Stores and gives access to four-momenta and invariant masses
/// \author Johannes Rauch, Daniel Greenwald
class FourMomenta : public StaticDataAccessor
{
public:

    /// Constructor
    FourMomenta(Model* m);

    /// check consistency
    bool consistent() const;

    /// Fill 4-momenta
    /// \param d DataPoint to fill
    /// \param dataPartitionIndex for status tracking
    virtual void calculate(DataPoint& d, unsigned dataPartitionIndex = 0) override;

    /// \name Getters
    /// @{

    /// Access 4-momenta (const)
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return 4-momentum of
    FourVector<double> p(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// Access invariant mass squared
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return squared mass of
    double m2(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
    { return pow(m(d, pc), 2); }

    /// Access invariant mass
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return mass of
    double m(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// \return initial-state four-momentum (const)
    /// \param d DataPoint to get data from
    const FourVector<double> initialStateMomentum(const DataPoint& d) const;

    /// \return vector of final-state four-momenta (const)
    /// \param d DataPoint to get data from
    const std::vector<FourVector<double> > finalStateMomenta(const DataPoint& d) const;

    /// \return masses
    std::shared_ptr<RealCachedDataValue> mass()
    { return M_; }

    /// \return masses (const)
    std::shared_ptr<RealCachedDataValue> mass() const
    { return M_; }

    /// \return momentum
    std::shared_ptr<FourVectorCachedDataValue> momentum()
    { return P_; }

    /// \return momentum (const)
    std::shared_ptr<FourVectorCachedDataValue> momentum() const
    { return P_; }

    /// @}

    /// print all masses
    std::ostream& printMasses(const DataPoint& d, std::ostream& os = std::cout) const;

    virtual std::string data_accessor_type() const override
    {return "FourMomenta"; }

    /// grant friend status to Model to call setFourMomenta and addParticleCombination
    friend class Model;

protected:

    /// set final-state four-momenta
    void setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, unsigned dataPartitionIndex = 0);

    /// looks for ISP when adding ParticleCombination's
    unsigned addParticleCombination(std::shared_ptr<ParticleCombination> pc) override;

    /// override to do nothing, since FourMomenta doesn't rely on parents being set.
    void pruneSymmetrizationIndices() override
    {}

private:

    /// Symmetrization index of initial state
    int ISPIndex_;

    /// Symmetrization indices of final states
    std::vector<int> FSPIndices_;

    /// four-vector of particle combinations
    std::shared_ptr<FourVectorCachedDataValue> P_;

    /// invariant mass of particle combinations [GeV]
    std::shared_ptr<RealCachedDataValue> M_;

};

}
#endif

