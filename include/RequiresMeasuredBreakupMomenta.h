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

#ifndef yap_RequiresMeasuredBreakupMomenta_h
#define yap_RequiresMeasuredBreakupMomenta_h

namespace yap {

/// \class RequiresMeasuredBreakupMomenta
/// \brief Base class to be inherited from to denote that an object
/// requires a model calculate measured breakup momenta
class RequiresMeasuredBreakupMomenta
{
public:

    /// Constructor
    /// \param r Whether to require measured breakup momenta
    RequiresMeasuredBreakupMomenta(bool r = true) : Requires_(r) {}

    /// \return Whether object requires measured breakup momenta be
    /// calculated by model
    bool requiresMeasuredBreakupMomenta() const
    { return Requires_; }

private:

    bool Requires_;

};

}
#endif
