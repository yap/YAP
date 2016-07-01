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
#include <vector>

namespace yap {

class DecayingParticle;
class FinalStateParticle;
class MassShape;
class Resonance;

/// \class lineInputIterator
/// Iterator to read an input stream line by line.
/// \attention The line is read when the iterator is increased.
/// \author Paolo Di Giglio
template <class StringT = std::string>
class LineInputIterator : public std::iterator<std::input_iterator_tag, StringT>
{
public:
	typedef typename StringT::value_type  char_type;
	typedef typename StringT::traits_type traits_type;
	typedef std::basic_istream<char_type, traits_type> istream_type;

	LineInputIterator(istream_type& is) : InputStream_(is) {};

	const StringT& operator*() const { return Value_; };
	const StringT* const operator->() const { return &Value_;}

	/// pre-increment operator (read line in)
	LineInputIterator<StringT>& operator++()
	{
		/// \todo Check for `getline()` return values
		getline(InputStream_, Value_);
	}

	/// post-increment operator (read line in)
	LineInputIterator<StringT> operator++(int)
	{
		LineInputIterator<StringT> prev(*this);
		// read line in
		++(*this);
		return prev;
	}
	
	friend const bool operator==(const LineInputIterator<StringT>& lhs, const LineInputIterator<StringT>& rhs)
	{ return (lhs.InputStream_ == rhs.InputStream_) && (lhs.Value_ == rhs.Value_); };
	
	friend const bool operator!=(const LineInputIterator<StringT>& lhs, const LineInputIterator<StringT>& rhs)
	{ return !(lhs == rhs); };

private:
	istream_type InputStream_;
	StringT      Value_;
};

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

    /// Constructor
    /// \param pdlFile Path to a pdl file like used by EvtGen
    ParticleFactory(const std::string pdlFile);

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


    /// add ParticleTableEntry to #particleTable_
    /// \param entry ParticleTableEntry to add to #particleTable_
    void addParticleTableEntry(ParticleTableEntry entry);

    // find PDG number by particle name
    // \return PDG code number
    // \param name Particle name as listed in particle table
    int pdgCode(std::string name) const;

    /// @}

private:
    /// read pdl file and fill #particleTable_
    void readPDT(const std::string fname);

    /// maps PDGCodes to ParticleTableEntry's
    std::map<int, ParticleTableEntry> particleTable_;

};

}

#endif
