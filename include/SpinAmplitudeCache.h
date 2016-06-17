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

#ifndef yap_SpinAmplitudeCache_h
#define yap_SpinAmplitudeCache_h

#include "fwd/Spin.h"

#include "Exceptions.h"
#include "SpinAmplitude.h"
#include "WeakPtrCache.h"

#include <memory>

namespace yap {

class Model;

/// \class SpinAmplitudeCache
/// \brief Caches SpinAmplitudes
/// \author Johannes Rauch, Daniel Greenwald
class SpinAmplitudeCache :
    public WeakPtrCache<SpinAmplitude>
{
public:

    /// Constructor
    /// \param model raw pointer to Model this cache belongs to
    SpinAmplitudeCache(Model* model = nullptr) :
        WeakPtrCache(), Model_(model) {}

    /// Equality
    bool equal(const std::shared_ptr<SpinAmplitude>& A, const std::shared_ptr<SpinAmplitude>& B) const override
    { return (A.get() == B.get()) or A->equalTo(*B); }

    /// retrieve or create SpinAmplitude
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector of daughters
    /// \param L orbital angular momentum
    /// \param two_S 2 * the total spin angular momentum
    std::shared_ptr<SpinAmplitude> spinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned L, unsigned two_S)
    {
        auto retVal = operator[](create(two_J, two_j, L, two_S));
        retVal->setModel(Model_);
        return retVal;
    }

    /// Check consistency of cache. Skips expired entries.
    bool consistent() const
    {
        bool C = true;
        for (auto it = begin(); it != end(); ++it)
            if (!it->expired())
                C &= it->lock()->consistent();
        return C;
    }

    /// grant friend status to Model to set itself owner
    friend class Model;

protected:

    /// set raw pointer to owning Model
    void setModel(Model* model)
    {
        if (Model_ != nullptr and Model_ != model)
            throw exceptions::Exception("Model already set.", "SpinAmplitudeCache::setModel");
        Model_ = model;
    }

private:

    /// override in inherting classes
    /// \return shared_ptr to SpinAmplitude object
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector of daughters
    /// \param L orbital angular momentum
    /// \param two_S 2 * the total spin angular momentum
    virtual std::shared_ptr<SpinAmplitude> create(unsigned two_J, const SpinVector& two_j, unsigned L, unsigned two_S) const = 0;

    /// raw pointer to Model this cache belongs to
    Model* Model_;

};

}

#endif
