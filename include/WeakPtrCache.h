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

#ifndef yap_WeakPtrCache_h
#define yap_WeakPtrCache_h

#include <memory>
#include <set>
#include <string>

namespace yap {

/// \class WeakPtrCache
/// \brief Template for a cache of weak_ptr's to objects
/// \author Daniel Greenwald

template <class T>
class WeakPtrCache
{
public:

    /// cache storage type
    using cache_type = std::set<std::weak_ptr<T>, std::owner_less<std::weak_ptr<T> > >;

    /// override to implement equivalence checking
    virtual bool equiv(const std::shared_ptr<T>& A, const std::shared_ptr<T>& B) const = 0;

    /// \name constructors and assignment operators
    /// @{

    /// Default constructor (defaulted)
    WeakPtrCache() = default;

    /// copy constructor (defaulted)
    WeakPtrCache(const WeakPtrCache&) = default;

    /// move constructor (defaulted)
    WeakPtrCache(WeakPtrCache<T>&&) = default;

    /// Construct cache from vector
    WeakPtrCache(std::vector<std::shared_ptr<T> > V)
    { for (auto& v : V) operator[](v); }

    /// virtual desctructor (defaulted)
    virtual ~WeakPtrCache() = default;

    /// copy assignment operator (defaulted)
    WeakPtrCache<T>& operator=(const WeakPtrCache<T>&) = default;

    /// move assignment operator (defaulted)
    WeakPtrCache<T>& operator=(WeakPtrCache<T>&&) = default;

    /// @}

    /// check if cache contains element equating to t
    /// \param t shared_ptr to object to search for equivalent of
    std::weak_ptr<T> find(std::shared_ptr<T> t) const
    {
        if (!t)
            return std::weak_ptr<T>();

        // search for equivalent
        auto it = std::find_if(Cache_.begin(), Cache_.end(), [&](const std::weak_ptr<T>& w) {return equiv(w.lock(), t);});

        if (it == Cache_.end())
            // if not found
            return std::weak_ptr<T>();

        return *it;
    }

    /// \return shared_ptr from Cache, if it exists, otherwise adds it to cache.
    /// \param t Shared ptr to object to retrieve from or add to cache
    std::shared_ptr<T> operator[](std::shared_ptr<T> t)
    {
        auto w = find(t);

        // if w is valid (t is found in cache)
        if (!w.expired())
            return w.lock();

        // else add to cache
        Cache_.emplace(t);
        return t;
    }

    /// remove expired Cache_ elements
    void removeExpired()
    {
        for (auto it = Cache_.begin(); it != Cache_.end(); ) {
            if (it->expired())
                it = Cache_.erase(it);
            else
                it++;
        }
    }

    /// \name access to cache
    /// @{

    /// \return iterator to begin
    typename cache_type::iterator begin()
    { return Cache_.begin(); }

    /// \return const_iterator to begin
    typename cache_type::const_iterator begin() const
    { return Cache_.begin(); }

    /// \return iterator to end
    typename cache_type::iterator end()
    { return Cache_.end(); }

    /// \return const_iterator to end
    typename cache_type::const_iterator end() const
    { return Cache_.end(); }

    /// @}

private:

    /// set of weak pointers to objects
    cache_type Cache_;

};

}

#endif
