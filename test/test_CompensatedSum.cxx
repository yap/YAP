#include <catch.hpp>

#include <CompensatedSum.h>

/**
 *  Test kahan summation algorithm
 */

TEST_CASE("CompensatedSum")
{
    // double has ~17 digits precision
    double init = 1.e17;

    yap::CompensatedSum<double> kahanSum(init);
    double normalSum(init);

    unsigned long N = 1000000;

    for (unsigned long i = 0; i < N; ++i) {
        kahanSum += 1;
        normalSum += 1;
    }

    // check that kahanSum gives the correct value
    REQUIRE( kahanSum == init + double(N) );

    // to make sure that kahanSum actually did something, require that the normal sum
    // did overflow and does NOT give the correct result
    REQUIRE( normalSum != init + double(N) );
}
