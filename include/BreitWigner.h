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

#ifndef yap_BreitWigner_h
#define yap_BreitWigner_h

#include "Amp.h"
#include "CalculationStatus.h"
#include "MassShape.h"

namespace yap {

class ParticleCombination;

/// \class BreitWigner
/// \brief Class for Breit-Wigner resonance shape
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is 1 / (mass^2 - s - i*mass*width)\n\n
/// Variables stored in #MassShapes::Parameters_:\n
///     MassShapes::#Parameters_[0] := nominal mass; set by setMass(double), returned by mass()\n
///     MassShapes::#Parameters_[1] := nominal width; set by setWidth(double), returned by width()\n\n
/// Values stored into DataPoint:
///     [0] := real(A(s))
///     [1] := imag(A(s))



class BreitWigner : public MassShape
{
public:

    /// \name Constructors, destructor, & operators
    /// @{

    /// Default constructor
    BreitWigner(InitialStateParticle* isp, double mass = -1, double width = -1);

    /// @}

    /// \name Getters
    /// @{

    /// \return nominal mass
    double mass() const
    { return Parameters_.at(0); }

    /// \return nominal squared masss
    double squaredMass() const
    { return mass() * mass(); }

    /// \return nominal width
    double width() const
    { return Parameters_.at(1); }

    /// @}

    /// \name Setters
    /// @{

    /// Set nominal mass
    void setMass(double m)
    { Parameters_.at(0) = m; }

    /// Set nominal width
    void setWidth(double w)
    { Parameters_.at(1) = w; }

    /// @}

    /// \name Amplitude related
    /// @{

    /// Calculate MassShape ampltude from squared mass;
    /// A = 1 / [Mass^2 - s - i * Mass * Width]
    /// \return amplitude evaluated at squared mass
    /// \param s squared mass to evaluate at
    virtual Amp calcAmplitudeS(double s);

    /// @}

    /// \name Bookkeeping related
    /// @{

    virtual bool consistent() const override;

    /// @}

protected:
    CalculationStatus CalcStatus_;

    Amp M2iMG_;                  // mass * mass - i * mass * width

};

}

#endif
