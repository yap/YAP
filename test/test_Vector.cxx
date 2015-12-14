#include <catch.hpp>

#include <Constants.h>
#include <FourVector.h>
#include <logging.h>
#include <NVector.h>
#include <ThreeVector.h>

#include <cmath>

TEST_CASE( "Vector" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);

    SECTION( "FourVector" ) {
        // cache J = 0, 1/2, 1 (out of order)
        REQUIRE( abs(yap::FourVector<double>({1,0,0,0})) == 1 );
        REQUIRE( abs(yap::FourVector<double>({1,1,0,0})) == 0 );
    }

}
