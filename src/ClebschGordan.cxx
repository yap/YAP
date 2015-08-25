#include "ClebschGordan.h"

#include "logging.h"

namespace yap {

//-------------------------
double clebschGordan(int two_j1, int two_m1, int two_j2, int two_m2, int two_J, int two_M)
{
    // check input parameters
    if (   not spinAndProjAreCompatible(two_j1, two_m1)
            or not spinAndProjAreCompatible(two_j2, two_m2)
            or not spinAndProjAreCompatible(two_J,  two_M )) {
        LOG(WARNING) << "spins and spin projections are inconsistent: "
                     << "(two_j1 = " << 0.5 * (two_j1) << ", two_m1 = " << 0.5 * (two_m1) << ", "
                     << "two_j2 = "  << 0.5 * (two_j2) << ", two_m2 = " << 0.5 * (two_m2) << ", "
                     << "two_J = "   << 0.5 * (two_J)  << ", two_M = "  << 0.5 * (two_M)  << ")" << std::endl;
        return 0;
    }
    if (not spinStatesCanCouple(two_j1, two_j2, two_J)) {
        LOG(WARNING) << "spins two_j1 = " << 0.5 * (two_j1) << " and two_j2 = " << 0.5 * (two_j2)
                     << " cannot couple to two_J = "  << 0.5 * (two_J) << std::endl;
        return 0;
    }
    if (two_m1 + two_m2 != two_M) {
        LOG(WARNING) << "spin projections two_m1 = " << 0.5 * (two_m1) << " and two_m2 = " << 0.5 * (two_m2)
                     << " cannot couple to two_M = " << 0.5 * (two_M) << std::endl;
        return 0;
    }

    double clebschVal = 0;
    // calculate function value and put intermediate values into cache
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

}
