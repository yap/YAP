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
/// Variables stored in #MassShapes::Parameters_:\n
///     MassShapes::#Parameters_[0] := nominal mass; set by setMass(double), returned by mass()\n
///     MassShapes::#Parameters_[1] := nominal width; set by setWidth(double), returned by width()\n\n
/// Amplitude is 1 / (mass^2 - s - i*mass*width)\n

class BreitWigner : public MassShape
{
public:

    /// \name Constructors, destructor, & operators
    /// @{

    /// Default constructor
    BreitWigner(double mass = -1, double width = -1);

    /// @}

    /// \name Getters
    /// @{

    /// \return nominal mass
    double mass() const
    { return Parameters_[0]; }

    /// \return mass-squared-dependent mass
    /// \param s squared mass to evaluate at
    double mass(double s) const
    { return mass(); }

    /// \return nominal squared masss
    double squaredmass() const
    { return mass() * mass(); }

    /// \return mass-squared-dependent squared mass
    /// \param s squared mass to evaluate at
    double squaredmass(double s) const
    { return squaredmass(); }

    /// \return nominal width
    double width() const
    { return Parameters_[1]; }

    /// \return mass-squared-dependent width
    /// \param s squared mass to evaluate at
    double width(double s) const
    { return width(); }

    /// @}

    /// \name Setters
    /// @{

    /// Set nominal mass
    void setMass(double m)
    { Parameters_[0] = m; }

    /// Set nominal width
    void setWidth(double w)
    { Parameters_[1] = w; }

    /// @}

    /// \name Amplitude related
    /// @{

    /// Calculate MassShape amplitude from DataPoint
    /// \return amplitude evaluated on DataPoint
    /// \param d DataPoint to evaluate on
    virtual Amp amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc) override;

    /// Calculate MassShape ampltude from squared mass;
    /// A = 1 / [Mass^2 - s - i * Mass * Width]
    /// \return amplitude evaluated at squared mass
    /// \param s squared mass to evaluate at
    virtual Amp amplitude(double s);

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
