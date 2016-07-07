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

#ifndef yap_DecayTreeVectorIntegralFwd_h
#define yap_DecayTreeVectorIntegralFwd_h

#include "fwd/DecayTree.h"

#include <complex>
#include <map>

namespace yap {

class DecayTreeVectorIntegral;

template <typename T>
struct IntegralElement;

/// \typedef RealIntegralElement
/// \ingroup Integration
using RealIntegralElement = IntegralElement<double>;

/// \typedef ComplexIntegralElement
/// \ingroup Integration
using ComplexIntegralElement = IntegralElement<std::complex<double> >;

/// \typedef DiagonalIntegralMap
/// maps shared_ptr to DecayTree to diagonal integral element
/// \ingroup Integration
using DiagonalIntegralMap = std::map<DecayTreeVector::value_type, RealIntegralElement>;

/// \typedef OffDiagonalIntegralMap
/// maps shared_ptr to DecayTree to off-diagonal integral element
/// \ingroup Integration
using OffDiagonalIntegralMap = std::map<std::array<DecayTreeVector::value_type, 2>, ComplexIntegralElement>;



}

#endif
