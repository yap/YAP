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

#ifndef yap_FreeAmplitude_h
#define yap_FreeAmplitude_h

#include "fwd/DecayChannel.h"
#include "fwd/SpinAmplitude.h"

#include "Constants.h"
#include "Parameter.h"

#include <complex>

namespace yap {

/// \class FreeAmplitude
/// \brief Stores complex free amplitude for the particular decay of a particle
/// \author Daniel Greenwald
/// \ingroup Parameters
///
/// A free amplitude is assigned for each decay of a particle into a
/// particular set of daughters (as indicated by the DecayChannel)
/// with a particular set of angular-momentum quantum numbers as given
/// by the SpinAmplitude and parent spin projection
class FreeAmplitude : public ComplexParameter
{
public:

    /// Constructor
    /// \param dc shared_ptr to DecayChannel
    /// \param sa shared_ptr to SpinAmplitude
    /// \param two_m (twice) the parent spin projection
    /// \param a value to initialize to
    FreeAmplitude(std::shared_ptr<DecayChannel> dc, std::shared_ptr<SpinAmplitude> sa,
                  int two_m, std::complex<double> a = Complex_1)
        : ComplexParameter(a), DecayChannel_(dc), SpinAmplitude_(sa), TwoM_(two_m) {}

    /// \return DecayChannel_
    const std::shared_ptr<DecayChannel>& decayChannel() const
    { return DecayChannel_; }

    /// \return SpinAmplitude_
    const std::shared_ptr<SpinAmplitude>& spinAmplitude() const
    { return SpinAmplitude_; }

    /// \return (twice) the spin projection of the decaying particle
    const int twoM() const
    { return TwoM_; }

private:

    /// DecayChannel for which this is a free amplitude
    std::shared_ptr<DecayChannel> DecayChannel_;

    /// SpinAmplitude for which this is a free amplitude
    std::shared_ptr<SpinAmplitude> SpinAmplitude_;

    /// (twice) spin projection of decaying particle
    int TwoM_;
};

}

#endif
