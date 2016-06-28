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

#ifndef yap_PHSP_h
#define yap_PHSP_h

#include "FourVector.h"
#include "Model.h"

#include <algorithm>
#include <random>
#include <vector>

namespace yap {

/// \return vector of four momenta for daughters uniformly randomly generated in phase space of model
template <class Generator>
const std::vector<FourVector<double> > phsp(const Model& M, const MassAxes& A, Generator& g, unsigned max_attempts = 1000)
{
    static std::uniform_real_distribution<double> uniform;

    // get mass^2 ranges:
    std::vector<std::array<double, 2> > r;
    r.reserve(A.size());
    std::transform(A.begin(), A.end(), std::back_inserter(r),
                   [&](const MassAxes::value_type& a){auto R = M.massRange(a, M.initialStateParticle()); R[0] *= R[0]; R[1] = R[1] * R[1] - R[0]; return R;});
    
    // create vector to store invariant masses in
    std::vector<double> m2(r.size(), 0);
    std::vector<FourVector<double> > P;

    unsigned n = 0;
    while (P.empty() and n < max_attempts) {
        // generate random point in hypercube of mass ranges
        std::transform(r.begin(), r.end(), m2.begin(), [&](const std::array<double, 2>& R) {return R[0] + R[1] * uniform(g);});
        P = M.calculateFourMomenta(A, m2, M.initialStateParticle());
    }
    return P;
}

}

#endif
