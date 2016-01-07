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

#ifndef yap_Exceptions_h
#define yap_Exceptions_h

#include <stdexcept>

namespace yap {

namespace exceptions {

/// \class Exception
/// \brief Base class for handling YAP exceptions
/// \defgroup Exceptions YAP exceptions
///
/// Use this class for creating exceptions with clarifying text, which
/// will not be caught by specific type. Use the other exception classes
/// inheriting from this, when an exception must be caught by specific
/// type.
struct Exception : public std::runtime_error {
public:
    Exception(const std::string& what_arg, const std::string& func_name)
        : std::runtime_error(what_arg + (func_name.empty() ? "" : " (" + func_name + ")")) {}
protected:
    Exception() : std::runtime_error("") {}
};

/// \class AngularMomentumNotConserved
/// \ingroup Exceptions
class AngularMomentumNotConserved : public Exception
{
};

/// \class InconsistentSpinProjection
/// \ingroup Exceptions
class InconsistentSpinProjection : public Exception {};

/// \class NonfiniteResult
/// \ingroup Exceptions
class NonfiniteResult : public Exception {};

}

}

#endif
