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

#ifndef  yap_pdl_h
#define  yap_pdl_h

#include "fwd/ParticleFactory.h"

#include "ParticleFactory.h"

#include <sstream>
#include <string>

namespace yap {

/// \class PDLIterator
/// Stream iterator targeted for `.pdl` files to read the input stream
/// line by line. It automatically discards comments and _set_-entries.
/// When either the EOF or the _end_ keyword is reached, the iterator
/// is set to its end state (i.e. `InputStream_` is set to `nullptr`).
/// \attention The line is read when the iterator is incremented.
/// \author Paolo Di Giglio
/// \ingroup ParticleFactory
class PDLIterator : public std::iterator<std::input_iterator_tag, std::string>
{
public:
    using char_type    = typename std::string::value_type;
    using traits_type  = typename std::string::traits_type;
    using istream_type = std::basic_istream<char_type, traits_type>;

    /// Construct and read the first line in
    PDLIterator(istream_type& is) : InputStream_(&is) { ++(*this); };

    /// Deference operator
    ///
    /// \return A pair of `int` (the PDG particle ID) and #ParticleTableEntry
    /// constructed from the file entry currently read in the `Value_` class
    /// member.
    /// \attention Isospin and parity are missing from `.pdl` format!
    const ParticleTableEntry operator*() const;

    /// Arrow iterator
    /// \return `i->m` is the same as `(*i).m`
    ParticleTableEntry operator->() const
    { return (*this).operator * (); };

    /// pre-increment operator (read line in)
    ///
    /// It ignores the comments (lines starting with '*') and
    /// all the lines starting with _set_. If the _end_ keyword
    /// is found or the EOF is reached, `InputStream_` is set
    /// to `nullptr`.
    PDLIterator& operator++();

    /// post-increment operator (read line in by calling `++(*this)`)
    PDLIterator operator++(int);

    /// Returns just an empty iterator, i.e. a default constructed one.
    static const PDLIterator& end()
    {
        static PDLIterator PDL_END;
        return PDL_END;
    }

    /// Check if pointers to streams are equal
    friend const bool operator==(const PDLIterator& lhs, const PDLIterator& rhs)
    { return lhs.InputStream_ == rhs.InputStream_; }

    /// Check if pointers to streams are not equal
    friend const bool operator!=(const PDLIterator& lhs, const PDLIterator& rhs)
    { return !(lhs == rhs); }

private:
    istream_type* InputStream_;
    std::string Value_;

    /// Default constructor (private).
    /// It's automatically called when the End Of File is reached
    PDLIterator() : InputStream_(nullptr) {};
};

/// Helper function to create a #ParticleFactory from an input `.pdl` file
/// \return A copy of the created #ParticleFactory
/// \ingroup ParticleFactory
ParticleFactory read_pdl_file(const std::string& filename);

}

#endif
