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
        // +=
        {
            yap::FourVector<double> v({1, 2, 3, 4});
            v += yap::FourVector<double>({5, 6, 7, 8});
            REQUIRE( v == yap::FourVector<double>({6, 8, 10, 12}) );
        }

        // -=
        {
            yap::FourVector<double> v({1, 2, 3, 4});
            v -= yap::FourVector<double>({5, 6, 7, 8});
            REQUIRE( v == yap::FourVector<double>({ -4, -4, -4, -4}) );
        }

        // *=
        {
            yap::FourVector<double> v({1, 2, 3, 4});
            v *= 2.;
            REQUIRE( v == yap::FourVector<double>({2, 4, 6, 8}) );
        }

        // +
        REQUIRE( yap::FourVector<double>({1, 2, 3, 4}) + yap::FourVector<double>({5, 6, 7, 8}) == yap::FourVector<double>({6, 8, 10, 12}) );

        // -
        REQUIRE( yap::FourVector<double>({1, 2, 3, 4}) - yap::FourVector<double>({5, 6, 7, 8}) == yap::FourVector<double>({ -4, -4, -4, -4}) );

        // *
        REQUIRE( yap::FourVector<double>({1, 2, 3, 4}) * yap::FourVector<double>({5, 6, 7, 8}) == -60 );

        // norm
        REQUIRE( norm(yap::FourVector<double>({1, 0, 0, 0})) == 1 );
        REQUIRE( norm(yap::FourVector<double>({5, 6, 7, 8})) == -124 );

        // abs
        REQUIRE( abs(yap::FourVector<double>({1, 0, 0, 0})) == 1 );
        REQUIRE( abs(yap::FourVector<double>({1, 1, 0, 0})) == 0 );

        // unit
        //REQUIRE( unit(yap::FourVector<double>({1,0,0,0})) == yap::FourVector<double>({1,0,0,0}) );

        // *

        // -

        // cross

        // angle

    }

}
