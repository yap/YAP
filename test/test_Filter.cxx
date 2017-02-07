#include <catch.hpp>

#include <Exceptions.h>
#include <Filter.h>

#include <vector>

TEST_CASE( "Filter" )
{

    std::vector<double> V(10, 1.);

    REQUIRE_THROWS_AS( yap::lone_elt(V), yap::exceptions::Exception );

}
