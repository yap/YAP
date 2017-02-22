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

#ifndef yap_set_parities_h
#define yap_set_parities_h

#include "ParticleTable.h"
#include "QuantumNumbers.h"

/// \function set_parities
/// set parities for selected mesons.
/// THIS IS NOT COMPLETE
void set_parities(ParticleTable& pdl) {

    // light unflavored mesons
    pdl["pi+"].quantumNumbers().setP(-1);
    pdl["pi-"].quantumNumbers().setP(-1);
    pdl["pi0"].quantumNumbers().setP(-1);

    pdl["eta"].quantumNumbers().setP(-1);

    pdl["f_0(500)"].quantumNumbers().setP(+1);

    pdl["rho+"].quantumNumbers().setP(-1);
    pdl["rho-"].quantumNumbers().setP(-1);
    pdl["rho0"].quantumNumbers().setP(-1);

    pdl["omega"].quantumNumbers().setP(-1);

    pdl["eta'"].quantumNumbers().setP(-1);

    pdl["f_0"].quantumNumbers().setP(+1);

    pdl["a_0+"].quantumNumbers().setP(+1);
    pdl["a_0-"].quantumNumbers().setP(+1);
    pdl["a_00"].quantumNumbers().setP(+1);

    pdl["phi"].quantumNumbers().setP(-1);

    pdl["h_1"].quantumNumbers().setP(+1);

    pdl["b_1+"].quantumNumbers().setP(+1);
    pdl["b_1-"].quantumNumbers().setP(+1);
    pdl["b_10"].quantumNumbers().setP(+1);

    pdl["a_1+"].quantumNumbers().setP(+1);
    pdl["a_1-"].quantumNumbers().setP(+1);
    pdl["a_10"].quantumNumbers().setP(+1);

    pdl["f_2"].quantumNumbers().setP(+1);

    pdl["f_1"].quantumNumbers().setP(+1);

    pdl["eta(2S)"].quantumNumbers().setP(-1);

    pdl["pi(2S)+"].quantumNumbers().setP(-1);
    pdl["pi(2S)-"].quantumNumbers().setP(-1);
    pdl["pi(2S)0"].quantumNumbers().setP(-1);

    pdl["a_2+"].quantumNumbers().setP(+1);
    pdl["a_2-"].quantumNumbers().setP(+1);
    pdl["a_20"].quantumNumbers().setP(+1);

    pdl["eta(1405)"].quantumNumbers().setP(-1);

    pdl["omega(2S)"].quantumNumbers().setP(-1);

    pdl["rho(2S)+"].quantumNumbers().setP(-1);
    pdl["rho(2S)-"].quantumNumbers().setP(-1);
    pdl["rho(2S)0"].quantumNumbers().setP(-1);

    pdl["eta(1475)"].quantumNumbers().setP(-1);

    pdl["f_0(1500)"].quantumNumbers().setP(+1);

    pdl["omega(1650)"].quantumNumbers().setP(-1);

    pdl["omega(3)(1670)"].quantumNumbers().setP(-1);


    // strange mesons
    pdl["K+"].quantumNumbers().setP(-1);
    pdl["K-"].quantumNumbers().setP(-1);

    pdl["K0"].quantumNumbers().setP(-1);
    pdl["anti-K0"].quantumNumbers().setP(-1);

    pdl["K_S0"].quantumNumbers().setP(-1);
    pdl["K_L0"].quantumNumbers().setP(-1);


    // charmed mesons
    pdl["D+"].quantumNumbers().setP(-1);
    pdl["D-"].quantumNumbers().setP(-1);

    pdl["D0"].quantumNumbers().setP(-1);
    pdl["anti-D0"].quantumNumbers().setP(-1);

}

#endif
