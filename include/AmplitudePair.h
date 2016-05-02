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

#ifndef yap_AmplitudePair_h
#define yap_AmplitudePair_h

#include "fwd/CachedDataValue.h"
#include "fwd/DecayChannel.h"
#include "fwd/Parameter.h"

#include "Constants.h"

#include <complex>
#include <memory>

namespace yap {

/// \name AmplitudePair
/// \brief Pair of a free an fixed amplitude of a DecayChannel
/// \author Johannes Rauch, Daniel Greenwald
struct AmplitudePair {
    AmplitudePair(DecayChannel* dc = nullptr, std::complex<double> free = Complex_1);
    std::shared_ptr<ComplexCachedDataValue> Fixed;
    std::shared_ptr<ComplexParameter> Free;
};

}

#endif
