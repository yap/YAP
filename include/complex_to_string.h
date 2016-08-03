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

#ifndef yap_complex_to_string_h
#define yap_complex_to_string_h

#include <complex>
#include <string>

namespace std {

/// \return string of complex number
template <typename T>
std::string to_string(const std::complex<T>& z)
{
    using std::to_string;
    return to_string(real(z))
           + (imag(z) >= 0 ? " + " + to_string(imag(z)) : " - " + to_string(imag(conj(z)))) + "i";
}

}

#endif
