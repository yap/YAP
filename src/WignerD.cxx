#include "WignerD.h"

#include "Constants.h"

namespace yap {

//-------------------------
double dFunctionCached::operator ()(const int j,
                                    const int m,
                                    const int n,
                                    const double theta)
{
    // check input parameters
    if (j >= (int)_maxJ) {
        LOG(ERROR) << "J = " << 0.5 * j << " is too large. maximum allowed J is "
                   << (_maxJ - 1) * 0.5 << ". Aborting..." << std::endl;
        throw;
    }
    if ((j < 0) or (abs(m) > j) or (abs(n) > j)) {
        LOG(ERROR) << "illegal argument for Wigner d^{J = " << 0.5 * j << "}"
                   << "_{M = " << 0.5 * m << ", M' = " << 0.5 * n << "}"
                   << "(theta = " << (theta) << "). Aborting..." << std::endl;
        throw;
    }

    // trivial case
    if (j == 0)
        return 1;

    // swap spin projections for negative angle
    int _m = m;
    int _n = n;
    double thetaHalf = theta / 2;
    if (theta < 0) {
        thetaHalf = abs(thetaHalf);
        std::swap(_m, _n);
    }

    const double cosThetaHalf = cos(thetaHalf);
    const double sinThetaHalf = sin(thetaHalf);

    double dFuncVal = 0;
    cacheEntryType*& cacheEntry = _cache[j][(j + _m) / 2][(j + _n) / 2];
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
        const int jpm = (j + _m) / 2;
        const int jpn = (j + _n) / 2;
        const int jmm = (j - _m) / 2;
        const int jmn = (j - _n) / 2;
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
            const int jmnk = j + (_m + _n) / 2 - 2 * k;
            const int jmk = (j + _m) / 2 - k;
            const int jnk = (j + _n) / 2 - k;
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
    for (int j = 0; j < (int)_maxJ; ++j)
        for (int m = -j; m <= j; m += 2)
            for (int n = -j; n <= j; n += 2) {
                cacheEntryType*& cacheEntry = _cache[j][(j + m) / 2][(j + n) / 2];
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
    for (int j = 0; j < (int)_maxJ; ++j)
        for (int m = -j; m <= j; m += 2)
            for (int n = -j; n <= j; n += 2)
                delete _cache[j][(j + m) / 2][(j + n) / 2];
}

//-------------------------
double dFunction(const int j,
                 const int m,
                 const int n,
                 const double theta) ///< Wigner d-function d^j_{m n}(theta)
{
    const double dFuncVal = dFunctionCached::instance()(j, m, n, theta);
    LOG(DEBUG) << "Wigner d^{J = " << 0.5 * j << "}_{M = " << 0.5 * m << ", "
               << "M' = " << 0.5 * n << "}(theta = " << (theta) << ") = "
               << (dFuncVal) << std::endl;
    return dFuncVal;
}

//-------------------------
Amp sphericalHarmonic(const int l,
                      const int m,
                      const double theta,
                      const double phi)
{
    // crude implementation using Wigner d-function
    const Amp YVal = sqrt((l + 1) / (4.*PI))
                     * std::exp(Amp(0, ((double)m / 2) * phi)) * dFunction(l, m, 0, theta);
    LOG(DEBUG) << "spherical harmonic Y_{l = " << 0.5 * l << "}^{m = " << 0.5 * m << "}"
               << "(phi = " << (phi) << ", theta = " << (theta) << ") = "
               << (YVal) << std::endl;
    return YVal;
}

//-------------------------
Amp DFunction(const int j,
              const int m,
              const int n,
              const double alpha,
              const double beta,
              const double gamma)
{
    const double arg = ((double)m / 2) * alpha + ((double)n / 2) * gamma;
    const Amp DFuncVal = std::exp(Amp(0, -arg)) * dFunction(j, m, n, beta);
    LOG(DEBUG) << "Wigner D^{J = " << 0.5 * j << "}_{M = " << 0.5 * m << ", "
               << "M' = " << 0.5 * n << "}(alpha = " << (alpha)
               << ", beta = " << (beta) << ", gamma = " << (gamma)
               << ") = " << (DFuncVal) << std::endl;
    return DFuncVal;
}

//-------------------------
Amp DFunctionConj(const int j,
                  const int m,
                  const int n,
                  const double alpha,
                  const double beta,
                  const double gamma)
{
    const Amp DFuncVal = std::conj(DFunction(j, m, n, alpha, beta, gamma));
    LOG(DEBUG) << "Wigner D^{J = " << 0.5 * j << "}*_{M = " << 0.5 * m << ", "
               << "M' = " << 0.5 * n << "}(alpha = " << (alpha)
               << ", beta = " << (beta) << ", gamma = " << (gamma)
               << ") = " << (DFuncVal) << std::endl;
    return DFuncVal;
}

}
