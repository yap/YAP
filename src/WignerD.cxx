#include "WignerD.h"

#include "Exceptions.h"
#include "logging.h"
#include "MathUtilities.h"
#include "spin.h"

#include <cmath>
#include <vector>

namespace yap {

namespace dMatrix {

// index is Kappa, in [0, min(J - M, J + N)]
using KappaFactorVector = std::vector<double>;

// first index is for J + M, in [0, twoJ]
// second index is for J + N, in [0, min(J + M, floor(J))]
using dMatrix = std::vector<std::vector<KappaFactorVector> >;

// Cache of d-matrix kappa term factors
// index is for (2J - 1), since J = 0 requires no cache
static std::vector<dMatrix> CachedMatrices_;

}

//-------------------------
const double dFunction(unsigned twoJ, int twoM, int twoN, double beta)
{
    // d^J_MN = (-)^(M-N) * d^J_NM
    if (twoN > twoM)
        return pow_negative_one((twoM - twoN) / 2) * dFunction(twoJ, twoN, twoM, beta);
    // N <= M now

    // d^J_MN = (-)^(M-N) * d^J_(-M)(-N)
    if (twoN > 0)
        return pow_negative_one((twoM - twoN) / 2) * dFunction(twoJ, -twoM, -twoN, beta);
    // N <= 0 now

    // check M
    if (is_odd(twoJ + twoM))
        throw exceptions::Exception(std::string("M = ") + spin_to_string(twoM)
                                    + " is inconsistent for J = " + spin_to_string(twoJ),
                                    "WignerD::dFunction");
    if (is_odd(twoJ + twoM))
        throw exceptions::Exception(std::string("N = ") + spin_to_string(twoN)
                                    + " is inconsistent for J = " + spin_to_string(twoJ),
                                    "WignerD::dFunction");

    // check that neither |M| nor |N| exceeds J
    // this returns 0 rather than throws because it is mathematically sound
    if (std::abs(twoM) > (int)twoJ or std::abs(twoN) > (int)twoJ)
        return 0;

    // trivial case of J = 0
    if (twoJ == 0)
        return 1;

    // cache dMatrix for J if necessary
    dMatrix::cache(twoJ);
    // if problem with caching (should not happen!)
    if (twoJ > dMatrix::CachedMatrices_.size())
        throw exceptions::Exception(std::string("could not cache Wigner d function for J = ") + spin_to_string(twoJ),
                                    "WignerD::dFunction");

    const dMatrix::KappaFactorVector& KF = dMatrix::CachedMatrices_[twoJ - 1][(twoJ + twoM) / 2][(twoJ + twoN) / 2];

    unsigned MminusN = (twoM - twoN) / 2;

    // sum over powers of cosine and sine of beta / 2 and multiply by factor
    double cosHalfBeta = cos(beta / 2); // power for cos is [2J - (M - N) - 2K]
    double sinHalfBeta = sin(beta / 2); // power for sin is [(M - N) + 2K]
    double dMatrixElement = 0;
    for (unsigned K = 0; K < KF.size(); ++K)
        dMatrixElement += KF[K] * pow(cosHalfBeta, twoJ - MminusN - 2 * K) * pow(sinHalfBeta, MminusN + 2 * K);

    return dMatrixElement;
}

//-------------------------
void dMatrix::cache(unsigned twoJ)
{
    /// d-matrix has already been cached for this spin
    if (twoJ == 0 or (twoJ <= CachedMatrices_.size() and !CachedMatrices_[twoJ - 1].empty()))
        return;

    /// resize d-matrix vector to hold up to spin J
    if (twoJ > CachedMatrices_.size())
        CachedMatrices_.resize(twoJ);

    double J = (double)twoJ / 2;

    dMatrix dJ(twoJ + 1);

    for (unsigned JplusM = 0; JplusM <= twoJ; ++JplusM) {

        unsigned JminusM = twoJ - JplusM;

        // = sqrt( (J + M)! * (J - M)! )
        double JMFactor = sqrt(std::tgamma(JplusM + 1) * std::tgamma(JminusM + 1));

        dJ[JplusM].resize(std::min(JplusM, (unsigned)std::floor(J)) + 1);

        for (unsigned JplusN = 0; JplusN < dJ[JplusM].size(); ++JplusN) {

            unsigned JminusN = twoJ - JplusN;
            unsigned MminusN = JplusM - JplusN;

            // = sqrt( (J + N)! * (J - N)! )
            double JNFactor = sqrt(std::tgamma(JplusN + 1) * std::tgamma(JminusN + 1));

            dJ[JplusM][JplusN].resize(std::min(JminusM, JplusN) + 1);

            // minK is 0, by choice that N <= M
            for (unsigned K = 0; K < dJ[JplusM][JplusN].size(); ++K)
                dJ[JplusM][JplusN][K] = pow_negative_one(K + MminusN) * JMFactor * JNFactor
                                        / std::tgamma(JminusM - K + 1) / std::tgamma(JplusN - K + 1) / std::tgamma(K + MminusN + 1) / std::tgamma(K + 1);
        }
    }
    CachedMatrices_[twoJ - 1] = dJ;
}

//-------------------------
const unsigned dMatrix::cacheSize()
{
    unsigned totSize = sizeof(CachedMatrices_);
    for (const auto& dJ : CachedMatrices_) {
        totSize += sizeof(dJ);
        for (const auto& dJrow : dJ) {
            totSize += sizeof(dJrow);
            for (const auto& dJelt : dJrow) {
                totSize += sizeof(dJelt) + dJelt.size() * sizeof(KappaFactorVector::value_type);
            }
        }
    }
    return totSize;
}

}
