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

#ifndef yap_BlattWeisskopfFwd_h
#define yap_BlattWeisskopfFwd_h

#include <map>
#include <memory>

namespace yap {

class BlattWeisskopf;

/// \typedef BlattWeisskopfMap
/// Map of unsigned (for angular momentum) to shared_ptr to BlattWeisskopf
using BlattWeisskopfMap = std::map<unsigned, std::shared_ptr<BlattWeisskopf> >;

}

#endif
