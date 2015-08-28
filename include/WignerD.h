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

#include "Amp.h"
#include "logging.h"
#include "MathUtilities.h"

#include <cmath>
#include <math.h>

namespace yap {

/// \class WignerD
/// \brief Wigner d and D functions
/// \author Daniel Greenwald, Johannes Rauch
/// This class has been copied from rootpwa and modified
//-------------------------------------------------------------------------
//
// Description:
// optimized Wigner d-function d^j_{m n}(theta) with caching
// used as basis for optimized spherical harmonics Y_l^m(theta, phi)
// as well as for optimized D-function D^j_{m n}(alpha, beta,
// gamma) and D-function in reflectivity basis
//
// based on PWA2000 function d_jmn_b() in pputil.cc
//
// !NOTE! spin j and projections m and n are in units of hbar/2
//
//
// Author List:
// Boris Grube TUM (original author)
//
//
//-------------------------------------------------------------------------

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


    ///< get singleton instance
    static dFunctionCached& instance() { return _instance; }

    ///< returns d^j_{m n}(theta)
    double operator ()(const int j,
                       const int m,
                       const int n,
                       const double theta);

    ///< returns caching flag
    static bool useCache()
    { return _useCache; }

    ///< sets caching flag
    static void setUseCache(const bool useCache = true)
    { _useCache = useCache; }

    ///< returns cache size in bytes
    static unsigned int cacheSize();



private:

    dFunctionCached () { }
    ~dFunctionCached();
    dFunctionCached (const dFunctionCached&);
    dFunctionCached& operator =(const dFunctionCached&);

    static dFunctionCached _instance; ///< singleton instance
    static bool _useCache; ///< if set to true, cache is used

    static const unsigned int _maxJ = 41; ///< maximum allowed angular momentum * 2 + 1
    static cacheEntryType* _cache[_maxJ][_maxJ + 1][_maxJ + 1]; ///< cache for intermediate terms [j][m][n]
};


dFunctionCached dFunctionCached::_instance;
bool dFunctionCached::_useCache = true;

dFunctionCached::cacheEntryType*
dFunctionCached::_cache[_maxJ][_maxJ + 1][_maxJ + 1];


///< Wigner d-function d^j_{m n}(theta)
double dFunction(const int j,
                 const int m,
                 const int n,
                 const double theta);


///< spherical harmonics Y_l^{m}(theta, phi)
Amp sphericalHarmonic
(const int l,
 const int m,
 const double theta,
 const double phi);

///< Wigner D-function D^j_{m n}(alpha, beta, gamma)
Amp DFunction(const int j,
              const int m,
              const int n,
              const double alpha,
              const double beta,
              const double gamma);


///< complex conjugate of Wigner D-function D^j_{m n}(alpha, beta, gamma)
Amp DFunctionConj(const int j,
                  const int m,
                  const int n,
                  const double alpha,
                  const double beta,
                  const double gamma);


}

#endif
