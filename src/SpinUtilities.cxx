#include "SpinUtilities.h"

#include "logging.h"
#include "MathUtilities.h"

#include <cmath>
#include <math.h>

namespace yap {

//-------------------------
double clebschGordan(int two_j1, int two_m1, int two_j2, int two_m2, int two_J, int two_M)
{
    // check input parameters
    if (   not spinAndProjAreCompatible(two_j1, two_m1)
            or not spinAndProjAreCompatible(two_j2, two_m2)
            or not spinAndProjAreCompatible(two_J,  two_M )) {
        LOG(DEBUG) << "spins and spin projections are inconsistent: "
                     << "(j1 = " << spinToString(two_j1) << ", m1 = " << spinToString(two_m1) << ", "
                     << "j2 = "  << spinToString(two_j2) << ", m2 = " << spinToString(two_m2) << ", "
                     << "J = "   << spinToString(two_J)  << ", M = "  << spinToString(two_M)  << ")" << std::endl;
        return 0;
    }
    if (not spinStatesCanCouple(two_j1, two_j2, two_J)) {
        LOG(DEBUG) << "spins j1 = " << spinToString(two_j1) << " and j2 = " << spinToString(two_j2)
                     << " cannot couple to J = "  << spinToString(two_J) << std::endl;
        return 0;
    }
    if (two_m1 + two_m2 != two_M) {
        LOG(DEBUG) << "spin projections m1 = " << spinToString(two_m1) << " and m2 = " << spinToString(two_m2)
                     << " cannot couple to M = " << spinToString(two_M) << std::endl;
        return 0;
    }

    double clebschVal = 0;
    // calculate function value
    int nu = 0;
    while (    ((two_j1 - two_j2 - two_M) / 2 + nu < 0)
               or ((two_j1 - two_m1)     / 2 + nu < 0))
        nu++;

    double sum = 0;
    int d1, d2, n1;
    while (     ((d1 = (two_J - two_j1 + two_j2) / 2 - nu) >= 0)
                and ((d2 = (two_J + two_M)       / 2 - nu) >= 0)
                and ((n1 = (two_j2 + two_J + two_m1) / 2 - nu) >= 0)) {
        const int d3 = (two_j1 - two_j2 - two_M) / 2 + nu;
        const int n2 = (two_j1 - two_m1)     / 2 + nu;
        sum +=   powMinusOne(nu + (two_j2 + two_m2) / 2) * factorial(n1) * factorial(n2)
                 / (  factorial(nu) * factorial(d1)
                      * factorial(d2) * factorial(d3));
        nu++;
    }

    if (sum == 0)
        return 0;

    const double N1 = factorial((two_J  + two_j1 - two_j2) / 2);
    const double N2 = factorial((two_J  - two_j1 + two_j2) / 2);
    const double N3 = factorial((two_j1 + two_j2 - two_J ) / 2);
    const double N4 = factorial((two_J + two_M) / 2);
    const double N5 = factorial((two_J - two_M) / 2);

    const double D0 = factorial((two_j1 + two_j2 + two_J) / 2 + 1);
    const double D1 = factorial((two_j1 - two_m1) / 2);
    const double D2 = factorial((two_j1 + two_m1) / 2);
    const double D3 = factorial((two_j2 - two_m2) / 2);
    const double D4 = factorial((two_j2 + two_m2) / 2);

    const double A  = (two_J + 1) * N1 * N2 * N3 * N4 * N5 / (D0 * D1 * D2 * D3 * D4);

    clebschVal = sqrt(A) * sum;

    return clebschVal;
}

//-------------------------
bool spinAndProjAreCompatible(const int spin, const int spinProj)
{
    return (spin >= 0) and (std::abs(spinProj) <= spin) and isEven(spin - spinProj);
}

//-------------------------
bool spinStatesCanCouple(const int two_j1,
                         const int two_j2,
                         const int two_J)
{
    if ((two_j1 < 0) or (two_j2 < 0))
        return false;

    // make sure J is in physical allowed range
    if (   (two_J < std::abs(two_j1 - two_j2))
            or (two_J > two_j1 + two_j2))
        return false;

    if (isOdd(two_j1 + two_j2 - two_J))  // check that J is in the half-integer or integer series, respectively
        return false;

    return true;
}

//-------------------------
std::string spinToString(const int twoJ)
{
    if (isEven(twoJ))
        return std::to_string(twoJ / 2);
    return std::to_string(twoJ) + "/2";
}


}
