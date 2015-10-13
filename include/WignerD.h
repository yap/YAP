/* YAP - Yet another PWA toolkit
 Copyright 2015, Technische Universitaet Muenchen,
 Authors: Daniel Greenwald, Johannes Rauch

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

#ifndef yap_WignerD_h
#define yap_WignerD_h

#include <complex>
#include <vector>

namespace yap {

/// Wigner d-function d^j_{two_m two_n}(theta)
double dFunction(const int two_j, const int two_m, const int two_n, const double theta);

/// spherical harmonics Y_l^{two_m}(theta, phi)
std::complex<double> sphericalHarmonic(const int two_l, const int two_m, const double theta, const double phi);

/// Wigner D-function D^j_{two_m two_n}(alpha, beta, gamma)
std::complex<double> DFunction(const int two_j, const int two_m, const int two_n, const double alpha, const double beta, const double gamma);

/// complex conjugate of Wigner D-function D^j_{two_m two_n}(alpha, beta, gamma)
std::complex<double> DFunctionConj(const int two_j, const int two_m, const int two_n, const double alpha, const double beta, const double gamma);


/// \class dFunctionCached
/// \brief cached Wigner d and D functions
/// \author Daniel Greenwald, Johannes Rauch
/// This class has been copied from rootpwa and modified
class dFunctionCached
{

public:

    struct cacheEntryType {
        cacheEntryType()
            : constTerm(0)
        { }

        double constTerm;
        std::vector<int> kmn1;
        std::vector<int> jmnk;
        std::vector<double> factor;
    };

    typedef std::vector<std::vector<std::vector<cacheEntryType> > > cacheType;


    /// get singleton instance
    static dFunctionCached& instance()
    { return _instance; }

    /// \return d^j_{two_m two_n}(theta)
    double operator ()(const int two_j, const int two_m, const int two_n, const double theta);

    /// \return caching flag
    static bool useCache()
    { return _useCache; }

    /// sets caching flag
    static void setUseCache(const bool useCache = true)
    { _useCache = useCache; }

    /// \return cache size in bytes
    static unsigned int cacheSize();


private:

    dFunctionCached () { }
    ~dFunctionCached();
    dFunctionCached (const dFunctionCached&);
    dFunctionCached& operator =(const dFunctionCached&);

    static dFunctionCached _instance; ///< singleton instance
    static bool _useCache; ///< if set to true, cache is used

    static const unsigned int _maxJ = 41; ///< maximum allowed angular momentum * 2 + 1
    static cacheEntryType* _cache[_maxJ][_maxJ + 1][_maxJ + 1]; ///< cache for intermediate terms [two_j][two_m][two_n]
};


}

#endif
