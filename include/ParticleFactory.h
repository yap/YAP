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

#ifndef yap_ParticleFactory_h
#define yap_ParticleFactory_h

#include "QuantumNumbers.h"

#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <vector>

namespace yap {

class DecayingParticle;
class FinalStateParticle;
class MassShape;
class Resonance;


/// \struct ParticleTableEntry
/// \brief Data container for storing particle information in database
/// \author Johannes Rauch, Daniel Greenwald
struct ParticleTableEntry : public QuantumNumbers {
    ParticleTableEntry(int pdg = 0, std::string name = "", QuantumNumbers q = QuantumNumbers(), double mass = -1, std::vector<double> parameters = {});
    bool consistent() const override;
    int PDG;
    std::string Name;
    double Mass;
    std::vector<double> MassShapeParameters;
};

/// \class ParticleFactory
/// \brief Factory class for easy creation of Particle objects from PDG codes.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle
class ParticleFactory
{
public:

    using ParticleTableMap = std::map<int, ParticleTableEntry>;
    using value_type = ParticleTableMap::value_type;
    using iterator   = ParticleTableMap::iterator;

    /// Constructor
    /// \param first First iterator to read the particle entries from
    /// \param last  Last iterator to read the particle entries from
    template <class Iter>
    ParticleFactory(Iter first, Iter last)
    { importTable(first, last); }

    /// Imports the table into the class from `first` to `last`
    /// \param first Iterator pointing to the first entry to import
    /// \param last  Iterator pointing to the last entry to import
    template <class Iter>
    void importTable(Iter first, Iter last)
    { std::copy(first, last, std::inserter(particleTable_, particleTable_.end())); }

    /// Create a FinalStateParticle from a PDG code
    /// \param PDG PDG code of particle to create
    /// \return shared pointer to new final state particle
    std::shared_ptr<FinalStateParticle> fsp(int PDG) const;

    /// Create an decayingParticle from a PDG code
    /// \param PDG PDG code of particle to create
    /// \param radialSize radial size of particle to create [GeV^-1]
    /// \return shared pointer to new DecayingParticle object
    std::shared_ptr<DecayingParticle> decayingParticle(int PDG, double radialSize) const;

    /// Create a Resonance from a PDG code and a MassShape
    /// \param PDG PDG code of particle to create
    /// \param radialSize Radial size of particle to create [GeV^-1]
    /// \param massShape Pointer to MassShape object describing resonance
    /// \return shared pointer to new Resonance object
    std::shared_ptr<Resonance> resonance(int PDG, double radialSize, std::shared_ptr<MassShape> massShape) const;

    /// Increment-assignment operator.
    ///
    /// Merges the #particleTable_ of the second operand to the
    /// #particleTable_ of the first operand.
    /// \param rhs The right hand side of the `+=`.
    /// \return A `const` reference to `*this` instance.
    const ParticleFactory& operator+=(const ParticleFactory& rhs);

    /// \name Particle table access
    /// @{

    /// get ParticleTableEntry from #particleTable_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    const ParticleTableEntry& particleTableEntry(int PDG) const;

    /// get ParticleTableEntry from #particleTable_ with safety checks
    /// \param name Name of particle in table
    const ParticleTableEntry& particleTableEntry(std::string name) const
    { return particleTableEntry(pdgCode(name)); }

    /// get #QuantumNumbers from #particleTable_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    const QuantumNumbers& quantumNumbers(int PDG) const
    { return static_cast<const QuantumNumbers&>(particleTableEntry(PDG)); }

    /// get #QuantumNumbers from #particleTable_ with safety checks
    /// \param name Name of particle in table
    const QuantumNumbers& quantumNumbers(std::string name) const
    { return static_cast<const QuantumNumbers&>(particleTableEntry(name)); }

    /// inserts the pair `int` and `ParticleTableEntry` to #particleTable_
    /// \param entry a A pair of PDG ID (integer) and ParticleTableEntry
    /// to add to #particleTable_
    std::pair<iterator, bool> insert(const value_type& entry);

    // find PDG number by particle name
    // \return PDG code number
    // \param name Particle name as listed in particle table
    int pdgCode(std::string name) const;

    /// @}

private:
    /// maps PDGCodes to ParticleTableEntry's
    std::map<int, ParticleTableEntry> particleTable_;
};

/// \class PDLIterator
/// Stream iterator targeted for `.pdl` files to read the input stream
/// line by line. It automatically discards comments and _set_-entries.
/// When either the EOF or the _end_ keyword is reached, the iterator
/// is set to its end state (i.e. `InputStream_` is set to `nullptr`).
/// \attention The line is read when the iterator is incremented.
/// \author Paolo Di Giglio
class PDLIterator : public std::iterator<std::input_iterator_tag, std::string>
{
public:
    using char_type    = typename std::string::value_type;
    using traits_type  = typename std::string::traits_type;
    using istream_type = std::basic_istream<char_type, traits_type>;

    /// Default constructor.
    /// It's automatically called when the End Of File is reached
    PDLIterator() : InputStream_(nullptr) {};

    /// Construct and read the first line in
    PDLIterator(istream_type& is) : InputStream_(&is) { ++(*this); };

    // XXX Needed?
//	const std::string& to_string() const { return Value_; };
//	const std::string* const operator->() const { return &Value_;}

    /// Deference operator
    ///
    /// \return A pair of `int` (the PDG particle ID) and #ParticleTableEntry
    /// constructed from the file entry currently read in the `Value_` class
    /// member.
    /// \attention Isospin and parity are missing from `.pdl` format!
    const ParticleFactory::value_type operator*() const
    {
        // create a stream with the entry string
        std::istringstream iss(Value_);

        // temporary variables
        std::string discard, particleName;
        int    stdhepid;
        double mass;
        double particleWidth;
        double particleMaxWidth; // not needed
        int    threeTimesCharge;
        int    twoTimesSpin;
        double lifetime; // dimension [c/mm]
        int    lundkc; // not needed

        // read the values in
        iss >> discard; // "add"
        iss >> discard; // "p"
        iss >> discard; // "Particle"
        iss >> particleName;
        iss >> stdhepid;
        iss >> mass;
        iss >> particleWidth;
        iss >> particleMaxWidth;
        iss >> threeTimesCharge;
        iss >> twoTimesSpin;
        iss >> lifetime;
        iss >> lundkc;

        return ParticleFactory::value_type(stdhepid, ParticleTableEntry(stdhepid, particleName,
                                           QuantumNumbers(twoTimesSpin, std::round(1. * threeTimesCharge / 3)),
                                           mass, {particleWidth}));
    }

    /// pre-increment operator (read line in)
    ///
    /// It ignores the comments (lines starting with '*') and
    /// all the lines starting with _set_. If the _end_ keyword
    /// is found or the EOF is reached, `InputStream_` is set
    /// to `nullptr`.
    PDLIterator& operator++()
    {
        // ignore comments and 'set' entries
        std::string command;
        do {
            // check if iterator reaches EOF
            if (InputStream_ && !getline(*InputStream_, Value_)) {
                InputStream_ = nullptr;
            } else {
                std::istringstream iss(Value_);
                iss >> command;
                // check if the 'end' keyword has been reached
                if (command == "end")
                    InputStream_ = nullptr;
            }
        } while (Value_[0] == '*' || command == "set");

        return *this;
    }

    /// post-increment operator (read line in by calling `++(*this)`)
    PDLIterator operator++(int)
    {
        PDLIterator prev(*this);
        // read line in
        ++(*this);
        return prev;
    }

    /// Returns just an empty iterator, i.e. a default constructed one.
    static PDLIterator end()
    {
        return PDLIterator();
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
};

ParticleFactory read_pdl_file(const std::string& filename);

}

#endif
