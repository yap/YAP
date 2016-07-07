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

#ifndef yap_VariableStatus_h
#define yap_VariableStatus_h

#include <string>

namespace yap {

/// \enum VariableStatus
enum class VariableStatus : int {
    changed   = -1,        ///< Parameter is free and has been changed
    fixed     = 0,         ///< Parameter is fixed
    unchanged = +1,        ///< Parameter is free but has not been changed
};

inline std::string to_string(const VariableStatus& s)
{
    switch (s) {
        case VariableStatus::changed:
            return "changed";
        case VariableStatus::fixed:
            return "fixed";
        case VariableStatus::unchanged:
            return "unchanged";
        default:
            return std::to_string(static_cast<int>(s));
    }
}

}

#endif
