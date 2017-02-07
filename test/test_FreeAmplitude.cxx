#include <catch.hpp>

#include <FreeAmplitude.h>

TEST_CASE( "FreeAmplitude" )
{

    REQUIRE_THROWS_AS( yap::FreeAmplitude(nullptr, nullptr, 1), yap::exceptions::Exception );

}
