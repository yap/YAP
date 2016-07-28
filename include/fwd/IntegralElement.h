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

#ifndef yap_IntegralElementFwd_h
#define yap_IntegralElementFwd_h

#include <complex>
#include <vector>

namespace yap {

template <typename T>
struct IntegralElement;

/// \typedef RealIntegralElement
/// \ingroup Integration
using RealIntegralElement = IntegralElement<double>;

/// \typedef RealIntegralElementVector
/// \ingroup Integration
using RealIntegralElementVector = std::vector<RealIntegralElement>;

/// \typedef ComplexIntegralElement
/// \ingroup Integration
using ComplexIntegralElement = IntegralElement<std::complex<double> >;

/// \typedef ComplexIntegralElementVector
/// \ingroup Integration
using ComplexIntegralElementVector = std::vector<ComplexIntegralElement>;

/// \typedef ComplexIntegralElementMatrix
/// \ingroup Integration
using ComplexIntegralElementMatrix = std::vector<ComplexIntegralElementVector>;

}

#endif
