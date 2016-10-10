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

#ifndef yap__AttributesFwd_h
#define yap__AttributesFwd_h

#include "fwd/AttributeUtilities.h"

namespace yap {

struct orbital_angular_momentum;
struct spin_angular_momentum;
struct spin_projection;

/// \typedef Functor class to check orbital angular momentum
using l_equals = check_attribute<orbital_angular_momentum>;

/// \typedef Functor class to check spin angular momentum
using s_equals = check_attribute<spin_angular_momentum>;
    
/// \typedef Functor class to check spin projection
using m_equals = check_attribute<spin_projection>;

struct to;
struct exactly_to;

struct is_fixed;
struct is_not_fixed;

struct has_free_amplitude;
struct has_decay_tree;
struct has_decay_channel;

struct parent_particle;
template <typename> class name_of;

struct has_a_mass;

/// \typedef functor to check particle name
using is_named = check_attribute<name_of<identity> >;

/// \typedef functor to return name of parent_particle
using parent_name = name_of<parent_particle>;

}

#endif
