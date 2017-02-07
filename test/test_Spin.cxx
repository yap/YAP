#include <catch.hpp>

#include <Spin.h>

TEST_CASE( "Spin" )
{

    auto SPV = yap::projections(2);
    REQUIRE_NOTHROW( yap::to_string(SPV) );

}
