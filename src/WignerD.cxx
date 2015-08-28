#include "WignerD.h"

#include "Constants.h"
#include "logging.h"
#include "MathUtilities.h"
#include "SpinUtilities.h"

#include <cmath>
#include <math.h>

namespace yap {

dFunctionCached dFunctionCached::_instance;
bool dFunctionCached::_useCache = true;

dFunctionCached::cacheEntryType*
dFunctionCached::_cache[_maxJ][_maxJ + 1][_maxJ + 1];


//-------------------------
double dFunctionCached::operator ()(const int two_j,
                                    const int two_m,
                                    const int two_n,
                                    const double theta)
{
    // check input parameters
    if (two_j >= (int)_maxJ) {
        LOG(ERROR) << "J = " << spinToString(two_j) << " is too large. maximum allowed J is "
                   << spinToString(_maxJ - 1) << ". Aborting..." << std::endl;
        throw;
    }
    if ((two_j < 0) or (abs(two_m) > two_j) or (abs(two_n) > two_j)) {
        LOG(ERROR) << "illegal argument for Wigner d^{J = " << spinToString(two_j) << "}"
                   << "_{M = " << spinToString(two_m) << ", M' = " << spinToString(two_n) << "}"
                   << "(theta = " << (theta) << "). Aborting..." << std::endl;
        throw;
    }

    // trivial case
    if (two_j == 0)
        return 1;

    // swap spin projections for negative angle
    int _m = two_m;
    int _n = two_n;
    double thetaHalf = theta / 2;
    if (theta < 0) {
        thetaHalf = abs(thetaHalf);
        std::swap(_m, _n);
    }

    const double cosThetaHalf = cos(thetaHalf);
    const double sinThetaHalf = sin(thetaHalf);

    double dFuncVal = 0;
    cacheEntryType*& cacheEntry = _cache[two_j][(two_j + _m) / 2][(two_j + _n) / 2];
    if (_useCache and cacheEntry) {
        // calculate function value using cache
        double sumTerm = 0;
        for (unsigned int i = 0; i < cacheEntry->factor.size(); ++i) {
            sumTerm += pow(cosThetaHalf, cacheEntry->kmn1[i])
                       * pow(sinThetaHalf, cacheEntry->jmnk[i]) / cacheEntry->factor[i];
        }
        dFuncVal = cacheEntry->constTerm * sumTerm;
    } else {
        // calculate function value and put intermediate values into cache
        if (_useCache)
            cacheEntry = new cacheEntryType();
        const int jpm = (two_j + _m) / 2;
        const int jpn = (two_j + _n) / 2;
        const int jmm = (two_j - _m) / 2;
        const int jmn = (two_j - _n) / 2;
        const double kk =   factorial(jpm)
                            * factorial(jmm)
                            * factorial(jpn)
                            * factorial(jmn);
        const double constTerm = powMinusOne(jpm) * sqrt(kk);
        if (cacheEntry)
            cacheEntry->constTerm = constTerm;

        double sumTerm = 0;
        const int mpn = (_m + _n) / 2;
        const int kMin = std::max(0, mpn);
        const int kMax = std::min(jpm, jpn);
        for (int k = kMin; k <= kMax; ++k) {
            const int kmn1 = 2 * k - (_m + _n) / 2;
            const int jmnk = two_j + (_m + _n) / 2 - 2 * k;
            const int jmk = (two_j + _m) / 2 - k;
            const int jnk = (two_j + _n) / 2 - k;
            const int kmn2 = k - (_m + _n) / 2;
            const double factor = ( factorial(k )
                                    * factorial(jmk )
                                    * factorial(jnk )
                                    * factorial(kmn2)) / powMinusOne(k);
            if (cacheEntry) {
                cacheEntry->kmn1.push_back (kmn1);
                cacheEntry->jmnk.push_back (jmnk);
                cacheEntry->factor.push_back(factor);
            }
            // using the 1 / factor here so that function value is the same as in PWA2000
            sumTerm += pow(cosThetaHalf, kmn1) * pow(sinThetaHalf, jmnk) / factor;
        }
        dFuncVal = constTerm * sumTerm;
    }

    return dFuncVal;
}

//-------------------------
unsigned int dFunctionCached::cacheSize() ///< returns cache size in bytes
{
    unsigned int size = _maxJ * (_maxJ + 1) * (_maxJ + 1) * sizeof(_cache[0][0][0]);
    for (int two_j = 0; two_j < (int)_maxJ; ++two_j)
        for (int two_m = -two_j; two_m <= two_j; two_m += 2)
            for (int two_n = -two_j; two_n <= two_j; two_n += 2) {
                cacheEntryType*& cacheEntry = _cache[two_j][(two_j + two_m) / 2][(two_j + two_n) / 2];
                if (cacheEntry) {
                    size += sizeof(*cacheEntry);
                    size += cacheEntry->kmn1.capacity () * sizeof(int);
                    size += cacheEntry->jmnk.capacity () * sizeof(int);
                    size += cacheEntry->factor.capacity() * sizeof(double);
                }
            }
    return size;
}

//-------------------------
dFunctionCached::~dFunctionCached()
{
    for (int two_j = 0; two_j < (int)_maxJ; ++two_j)
        for (int two_m = -two_j; two_m <= two_j; two_m += 2)
            for (int two_n = -two_j; two_n <= two_j; two_n += 2)
                delete _cache[two_j][(two_j + two_m) / 2][(two_j + two_n) / 2];
}

//-------------------------
double dFunction(const int two_j,
                 const int two_m,
                 const int two_n,
                 const double theta) ///< Wigner d-function d^j_{two_m two_n}(theta)
{
    const double dFuncVal = dFunctionCached::instance()(two_j, two_m, two_n, theta);
    LOG(DEBUG) << "Wigner d^{J = " << spinToString(two_j) << "}_{M = " << spinToString(two_m) << ", "
               << "M' = " << spinToString(two_n) << "}(theta = " << (theta) << ") = "
               << (dFuncVal) << std::endl;
    return dFuncVal;
}

//-------------------------
Amp sphericalHarmonic(const int two_l,
                      const int two_m,
                      const double theta,
                      const double phi)
{
    // crude implementation using Wigner d-function
    const Amp YVal = sqrt((two_l + 1) / (4.*PI))
                     * std::exp(Amp(0, ((double)two_m / 2) * phi)) * dFunction(two_l, two_m, 0, theta);
    LOG(DEBUG) << "spherical harmonic Y_{l = " << spinToString(two_l) << "}^{m = " << spinToString(two_m) << "}"
               << "(phi = " << (phi) << ", theta = " << (theta) << ") = "
               << (YVal) << std::endl;
    return YVal;
}

//-------------------------
Amp DFunction(const int two_j,
              const int two_m,
              const int two_n,
              const double alpha,
              const double beta,
              const double gamma)
{
    const double arg = ((double)two_m / 2) * alpha + ((double)two_n / 2) * gamma;
    const Amp DFuncVal = std::exp(Amp(0, -arg)) * dFunction(two_j, two_m, two_n, beta);
    LOG(DEBUG) << "Wigner D^{J = " << spinToString(two_j) << "}_{M = " << spinToString(two_m) << ", "
               << "M' = " << spinToString(two_n) << "}(alpha = " << (alpha)
               << ", beta = " << (beta) << ", gamma = " << (gamma)
               << ") = " << (DFuncVal) << std::endl;
    return DFuncVal;
}

//-------------------------
Amp DFunctionConj(const int two_j,
                  const int two_m,
                  const int two_n,
                  const double alpha,
                  const double beta,
                  const double gamma)
{
    const Amp DFuncVal = std::conj(DFunction(two_j, two_m, two_n, alpha, beta, gamma));
    LOG(DEBUG) << "Wigner D^{J = " << spinToString(two_j) << "}*_{M = " << spinToString(two_m) << ", "
               << "M' = " << spinToString(two_n) << "}(alpha = " << (alpha)
               << ", beta = " << (beta) << ", gamma = " << (gamma)
               << ") = " << (DFuncVal) << std::endl;
    return DFuncVal;
}

}
