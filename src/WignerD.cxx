#include "WignerD.h"

#include "logging.h"
#include "MathUtilities.h"
#include "SpinUtilities.h"

INITIALIZE_EASYLOGGINGPP

#include <cmath>
#include <math.h>
#include <vector>

namespace yap {

namespace dMatrix {

// index is Kappa, in [0, min(J - M, J + N)]
using KappaFactorVector = std::vector<double>;

// first index is for J + M, in [0, twoJ]
// second index is for J + N, in [0, min(J + M, floor(J))]
using dMatrix = std::vector<std::vector<KappaFactorVector> >;

// Cache of d-matrix kappa term factors
// index is for twoJ = twice the spin of the representation
static std::vector<dMatrix> CachedMatrices_;

// // Fill cache at initialization of code
// LOG(INFO) << "filling d-Matrix cache ...";
// for (unsigned char twoJ = 0; twoJ <= 4; ++twoJ)
//     cacheMatrix(twoJ);
// LOG(INFO) << "... done filling d-Matrix cache, size = " << cacheSize() << " for " << 4 << " cached matrices.";

}

//-------------------------
std::complex<double> DFunction(unsigned char twoJ, char twoM, char twoN, double alpha, double beta, double gamma)
{
    return std::exp(-Complex_i * (alpha * twoM + gamma * twoN) / 2.) * dFunction(twoJ, twoM, twoN, beta);
}

//-------------------------
double dFunction(unsigned char twoJ, char twoM, char twoN, double beta)
{
    // d^J_MN(-beta) = dJ_NM(beta)
    if (beta < 0)
        return dFunction(twoJ, twoN, twoM, -beta);
    // beta is now positive

    // d^J_MN = (-)^(M-N) * d^J_NM
    if (twoN > twoM)
        return powMinusOne((twoM - twoN) / 2) * dFunction(twoJ, twoN, twoM, beta);
    // N <= M now

    // d^J_MN = (-)^(M-N) * d^J_(-M)(-N)
    if (twoN > 0)
        return powMinusOne((twoM - twoN) / 2) * dFunction(twoJ, -twoM, -twoN, beta);
    // N <= 0 now

    // check M
    if (isOdd(twoM) != isOdd(twoJ)) {
        FLOG(ERROR) << "helicity M = " << spinToString(twoM) << " invalid for spin J = " << spinToString(twoJ);
        throw std::invalid_argument("twoM");
    }
    if (std::abs(twoM) > twoJ) {
        FLOG(WARNING) << "helicity M = " << spinToString(twoM) << " is larger than spin J = " << spinToString(twoJ) << "; matrix element is zero";
        return 0;
    }

    // check N
    if (isOdd(twoN) != isOdd(twoJ)) {
        FLOG(ERROR) << "helicity N = " << spinToString(twoN) << " invalid for spin J = " << spinToString(twoJ);
        throw std::invalid_argument("twoN");
    }
    if (std::abs(twoN) > twoJ) {
        FLOG(WARNING) << "helicity N = " << spinToString(twoN) << " is larger than spin J = " << spinToString(twoJ) << "; matrix element is zero";
        return 0;
    }

    // trivial case of J = 0
    if (twoJ == 0)
        return 1;

    // cache dMatrix for J if necessary
    dMatrix::cache(twoJ);
    // if problem with caching (should not happen!)
    if (twoJ >= dMatrix::CachedMatrices_.size()) {
        FLOG(ERROR) << "d matrix could not be cached for spin J = " << spinToString(twoJ);
        FLOG(ERROR) << "CachedMatrices_.size() = " << dMatrix::CachedMatrices_.size();
        throw;
    }

    const dMatrix::KappaFactorVector& KF = dMatrix::CachedMatrices_[twoJ][(twoJ + twoM) / 2][(twoJ + twoN) / 2];

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
void dMatrix::cache(unsigned char twoJ)
{
    /// d-matrix has already been cached for this spin
    if (twoJ < CachedMatrices_.size() or CachedMatrices_.empty() or !CachedMatrices_[twoJ].empty())
        return;

    /// resize d-matrix vector to hold up to spin J
    if (twoJ >= CachedMatrices_.size())
        CachedMatrices_.resize(twoJ + 1);

    double J = (double)twoJ / 2;

    dMatrix dJ((unsigned)std::floor(J) + 1);

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

            dJ[JplusM][JplusN].resize(std::min(JminusM, JplusN));

            // minK is 0, by choice that N <= M
            for (unsigned K = 0; K <= dJ[JplusM][JplusN].size(); ++K)
                dJ[JplusM][JplusN][K] = powMinusOne(K + MminusN) * JMFactor * JNFactor
                                        / std::tgamma(JminusM - K + 1) / std::tgamma(JplusN - K + 1) / std::tgamma(K + MminusN) / std::tgamma(K + 1);
        }
    }
    CachedMatrices_[twoJ] = dJ;
}

//-------------------------
unsigned int dMatrix::cacheSize()
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
