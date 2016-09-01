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

#ifndef yap_TwoBodyPhaseSpaceFactor_h
#define yap_TwoBodyPhaseSpaceFactor_h

#include "fwd/TwoBodyPhaseSpaceFactor.h"

#include "fwd/DataPoint.h"
#include "fwd/DecayChannel.h"
#include "fwd/MassShape.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/SpinAmplitude.h"

#include "PhaseSpaceFactor.h"
#include "RequiresMeasuredBreakupMomenta.h"

namespace yap {

/// Calculates a two-body phase-space factor: rho := 2 q / m
/// \ingroup PhaseSpaceFactor
/// \author Daniel Greenwald
class TwoBodyPhaseSpaceFactor :
    public PhaseSpaceFactor,
    public RequiresMeasuredBreakupMomenta
{
public:
    /// Constructor
    TwoBodyPhaseSpaceFactor(const Model& m);

    /// \return MeasuredBreakupMomenta::rho(d, pc)
    virtual const std::complex<double> value(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// get raw pointer to Model (const)
    virtual const Model* model() const
    { return Model_; }

    /// do nothing
    void calculate(DataPartition& D) const {}

    /// do nothing
    virtual void updateCalculationStatus(StatusManager& D) const  {}

protected:

    /// add ParticleCombination to SymmetrizationIndices_
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> pc);

private:
    
    /// raw pointer to owning model
    const Model* Model_;

};

/// Factory for creating TwoBodyPhaseSpaceFactor objects
/// \ingroup PhaseSpaceFactor  
/// \author Daniel Greenwald
class TwoBodyPhaseSpaceFactorFactory : public PhaseSpaceFactorFactory
{
private:
    /// private constructor; singleton
    TwoBodyPhaseSpaceFactorFactory() = default;

public:

    /// \return singleton instance
    static std::shared_ptr<PhaseSpaceFactorFactory> instance()
    {
        static auto Instance_ = std::shared_ptr<PhaseSpaceFactorFactory>(new TwoBodyPhaseSpaceFactorFactory());
        return Instance_;
    }

protected:

    /// \return TwoBodyPhaseSpaceFactor for two-body decay; nullptr otherwise
    virtual std::shared_ptr<PhaseSpaceFactor> phaseSpaceFactor(const DecayChannel& dc, const SpinAmplitude& sa, std::shared_ptr<MassShape> ms);

private:

    /// single instance of TwoBodyPhaseSpaceFactor to be used by all objects
    std::map<const Model*, std::shared_ptr<TwoBodyPhaseSpaceFactor> > PHSPFactors_;

};
  
}

#endif
