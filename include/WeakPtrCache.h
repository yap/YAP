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

#include <algorithm>
#include <memory>
#include <ostream>
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

    /// \name helper types
    /// @{

    /// object type
    using type = T;

    /// non-const object type
    using non_const_type = typename std::remove_const<type>::type;

    /// shared_ptr_type
    /// \brief std::shared_ptr to T
    using shared_ptr_type = std::shared_ptr<T>;

    /// weak_ptr_type
    /// \brief std::weak_ptr to T
    using weak_ptr_type = std::weak_ptr<T>;

    /// cache storage type
    /// \brief A std::set of weak_ptr_type
    using cache_type = std::set<weak_ptr_type, std::owner_less<weak_ptr_type> >;

    /// @}

    /// override to implement equality checking
    virtual bool equal(const shared_ptr_type& A, const shared_ptr_type& B) const = 0;

    /// \name constructors and assignment operators
    /// @{

    /// Default constructor (defaulted)
    WeakPtrCache() = default;

    /// copy constructor (defaulted)
    WeakPtrCache(const WeakPtrCache&) = default;

    /// move constructor (defaulted)
    WeakPtrCache(WeakPtrCache&&) = default;

    /// Construct cache from vector
    WeakPtrCache(std::vector<shared_ptr_type> V)
    { for (auto& v : V) operator[](v); }

    /// virtual desctructor (defaulted)
    virtual ~WeakPtrCache() = default;

    /// copy assignment operator (defaulted)
    WeakPtrCache& operator=(const WeakPtrCache&) = default;

    /// move assignment operator (defaulted)
    WeakPtrCache& operator=(WeakPtrCache&&) = default;

    /// @}

    /// check if cache contains element equating to t
    /// \param t shared_ptr to object to search for equivalent of
    weak_ptr_type find(shared_ptr_type t) const
    {
        if (!t)
            return weak_ptr_type();

        // search for equivalent
        auto it = std::find_if(Cache_.begin(), Cache_.end(), [&](const weak_ptr_type & w) {return !w.expired() and equal(w.lock(), t);});

        if (it == Cache_.end())
            // if not found
            return weak_ptr_type();

        return *it;
    }

    /// \return shared_ptr from Cache, if it exists, otherwise adds it to cache.
    /// \param t Shared ptr to object to retrieve from or add to cache
    shared_ptr_type operator[](shared_ptr_type t)
    {
        auto w = find(t);

        // if w is valid (t is found in cache)
        if (!w.expired())
            return w.lock();

        // else add to cache
        addToCache(t);
        return t;
    }

    /// \return whether empty
    bool empty() const
    { return Cache_.empty(); }

    /// \return size
    size_t size() const
    { return Cache_.size(); }

    /// \return number of expired Cache_ elements
    size_t count_expired() const
    { return std::count_if(Cache_.begin(), Cache_.end(), [](const weak_ptr_type & w) {return w.expired();}); }

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

protected:

    /// add element to cache
    virtual void addToCache(shared_ptr_type t)
    { Cache_.emplace(t); }

private:

    /// set of weak pointers to objects
    cache_type Cache_;

};

/// streamer
template <class T>
inline std::string to_string(const WeakPtrCache<T>& C)
{
    std::string s = "contains " + std::to_string(C.size()) + " elements, of which "
        + std::to_string(C.count_expired()) + " have expired";
    using std::to_string;
    for (const auto& w : C)
        if (!w.expired())
            s += "\n" + to_string(*w.lock());
    return s;
}

}

#endif
