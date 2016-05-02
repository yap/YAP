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

#ifndef yap_DecayTree_h
#define yap_DecayTree_h

#include "fwd/DecayTree.h"

#include "AmplitudePair.h"
#include "CachedDataValue.h"
#include "DecayChannel.h"

#include <algorithm>
#include <string>
#include <vector>

namespace yap {

/// \class DecayTree
/// \brief Class holding vectors of fixed and free amplitudes that define a decay tree
/// \author Johannes Rauch, Daniel Greenwald
class DecayTree
{
public:

    /// default constructor
    DecayTree() {}

    /// constructor
    DecayTree(const std::vector<AmplitudePair>& amplitudes) :
        Amplitudes_(amplitudes)
    {}

    /// constructor
    DecayTree(const AmplitudePair& ap) :
        Amplitudes_(1, ap)
    {}

    /// constructor
    DecayTree(const AmplitudePair& ap, const DecayTree& tree) :
        Amplitudes_()
    {
        Amplitudes_.reserve(tree.Amplitudes_.size() + 1);
        Amplitudes_.push_back(ap);
        Amplitudes_.insert(Amplitudes_.end(), tree.Amplitudes_.begin(), tree.Amplitudes_.end());
    }

    /// \name Getters
    /// @{

    const std::vector<AmplitudePair>&  amplitudes() const
    { return Amplitudes_; }

    std::vector<AmplitudePair>&  amplitudes()
    { return Amplitudes_; }

    /// @}


private:

    std::vector<AmplitudePair> Amplitudes_;

};

/// convert to string
inline std::string to_string(const DecayTree& t)
{
    std::string s("DecayTree with DecayChannels ");
    for (auto& ap : t.amplitudes()) {
        if (dynamic_cast<DecayChannel*>(ap.Fixed->owner()))
            s += to_string(*static_cast<DecayChannel*>(ap.Fixed->owner()));
        else
            s += "???";
        s += "; ";
    }
    s.erase(s.size() - 2, 2);

    return s;
}

}

#endif
