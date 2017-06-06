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

#include "fwd/FinalStateParticle.h"

#ifndef yap_DalitzPhspIntegral_h
#define yap_DalitzPhspIntegral_h

namespace yap {

double dalitz_phasespace_volume(double isp_mass, const FinalStateParticleVector& fsps, const double relErr = 1.e-7);

class DalitzIntegrand {
public:
    DalitzIntegrand(double isp_mass, const FinalStateParticleVector& fsps);

    double operator()(double mab2) const;

private:
    double M2_;
    double ma2_;
    double mb2_;
    double mc2_;
};

}

#endif
