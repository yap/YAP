#include "ClebschGordan.h"

#include <algorithm>
#include <cmath>

namespace yap {

//-------------------------
double ClebschGordan::coefficient(unsigned two_j1, int two_m1, unsigned two_j2, int two_m2, unsigned two_J, int two_M)
{
    // check input spin-projection compatibilities
    if (!consistent(two_j1, two_m1) or !consistent(two_j2, two_m2) or !consistent(two_J,  two_M ))
        return 0;

    // check whether input spins couple
    if (!couple(two_j1, two_j2, two_J))
        return 0;

    // check input spin-projection coupling
    if (two_M != two_m1 + two_m2)
        return 0;

    // z range dictated by factorials in denominator ( 1/n! = 0 when n < 0)
    unsigned z_min = std::max({0, two_j2 - two_m1 - J, two_j1 + two_m2 - two_J}) / 2;
    unsigned z_max = std::min({two_j1 + two_j2 - two_J, two_j1 - two_m1, two_j2 + two_m2}) / 2;
    
    double z_sum = 0;
    for (unsigned z = z_min; z <= z_max; ++z) {
        // z'th term := (-)^z / z! / (j1+j2-J-z)! / (j1-m1-z)! / (j2+m2-z)! / (J-j2+m1+z)! / (J-j1-m2+z)!
        z_sum += pow_negative_one(z)
            / std::tgamma(1 + z)
            / std::tgamma(1 + (two_j1 + two_j2 - two_J) / 2 - z)
            / std::tgamma(1 + (two_j1 - two_m1) / 2 - z)
            / std::tgamma(1 + (two_j2 + two_m2) / 2 - z)
            / std::tgamma(1 + (two_J - two_j2 + two_m1) / 2 + z)
            / std::tgamma(1 + (two_J - two_j1 - two_m2) / 2 + z);
    }

    // C-G coef = sqrt( (2J+1)! (j1+j2-J)! (j1-j2+J)! (j2-j1+J)! / (J+j1+j2+1)! )
    //          * sqrt( (j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)! (J+M)! (J-M)! )
    //          * z_sum
    return z_sum * sqrt(std::tgamma(1 + two_J + 1)
                        * std::tgamma(1 + (two_j1 + two_j2 - two_J) / 2)
                        * std::tgamma(1 + (two_j1 - two_j2 + two_J) / 2)
                        * std::tgamma(1 + (two_j2 - two_j1 + two_J) / 2)
                        / std::tgamma(1 + (two_j1 + two_j2 + two_J) / 2 + 1)
                        * std::tgamma(1 + (two_j1 + two_m1) / 2)
                        * std::tgamma(1 + (two_j1 - two_m1) / 2)
                        * std::tgamma(1 + (two_j2 + two_m2) / 2)
                        * std::tgamma(1 + (two_j2 - two_m2) / 2)
                        * std::tgamma(1 + (two_J + two_M) / 2)
                        * std::tgamma(1 + (two_J - two_M) / 2));
}

}
