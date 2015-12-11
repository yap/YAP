/* YAP - Yet another PWA toolkit
 Copyright 2015, Technische Universitaet Muenchen,
 Authors: Daniel Greenwald, Johannes Rauch

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/// \file
/// \brief Functions for caching and calculating Wigner d and D functions
/// \author Daniel Greenwald, Johannes Rauch
///
/// Calculation is based on M. E. Rose's _Elementary Theory of Angular Momentum_ (1957):
///
/// \f$ D^{J}_{MN}(\alpha, \beta, \gamma) \equiv \exp(-iM\alpha) d^{J}_{MN}(\beta) \exp(-iN\gamma)\f$
///
/// \f$ d^{J}_{MN}(\beta) \equiv \sum_{K} A (-)^{K+M-N} B_{K} (\cos\frac{\beta}{2})^{2J+N-M-2K} (\sin\frac{\beta}{2})^{M-N+2K} \f$
///
/// with \f$ A \equiv (J+N)!(J-N)!(J+M)!(J-M)!\f$
/// and \f$ B_{K} \equiv (J-M-K)!((J+N-K)!(K+M-N)!K!\f$.
///
/// The limits on K are governed by the factorials,
/// since \f$ (1 / X!) = 0 \f$ if \f$ X < 0\f$.
///
/// We only cache matrix elements with M in [-J, J] and N in [-J, min(0, m)].
/// This amounts to the lower triangle and the diagonal without the bottom right corner.
/// The uncached matrix elements are given by the by the symmetries
///   - \f$ d^{J}_{MN}(\beta) = (-)^(M-N) d^{J}_{NM}(\beta)\f$
///   - \f$ d^{J}_{MN}(\beta) = (-)^{M-N) d^{J}_{-N-M}(\beta)\f$

#ifndef yap_WignerD_h
#define yap_WignerD_h

#include "Constants.h"

#include <complex>

namespace yap {

/// Wigner D-function \f$ D^{J}_{M N}(\alpha, \beta, \gamma) \f$
std::complex<double> DFunction(unsigned char twoJ, char twoM, char twoN, double alpha, double beta, double gamma);

/// \return Wigner d-function \f$ d^{J}_{M N}(\beta) \f$
/// \param twoJ twice the total spin of system
/// \param twoM twice the first spin projection
/// \param twoN twice the second spin projection
/// \param beta rotation angle
double dFunction(unsigned char twoJ, char twoM, char twoN, double beta);


namespace dMatrix {

/// Cache d-matrix elements for representation of spin J
/// \param twoJ twice the spin of the representation
void cache(unsigned char twoJ);

/// \return cache size in bytes
unsigned cacheSize();

}

}

#endif
