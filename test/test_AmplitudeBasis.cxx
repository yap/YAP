#include <catch.hpp>

#include <AmplitudeBasis.h>

TEST_CASE( "AmplitudeBasis" )
{

    yap::amplitude_basis::canonical<double> C(1, 2, 3);
    REQUIRE_THROWS_AS( C[3], yap::exceptions::Exception );

}
