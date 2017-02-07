#include <catch.hpp>

#include <MathUtilities.h>

TEST_CASE( "MathUtilities" ) {

    REQUIRE(yap::rad(180.) == Approx(yap::pi()));

}
