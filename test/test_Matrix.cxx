#include <catch.hpp>

#include <Constants.h>
#include <Matrix.h>
#include <logging.h>

#include <cmath>

TEST_CASE( "Matrix" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);

    yap::ThreeMatrix<double> zero({0,0,0, 0,0,0, 0,0,0});
    yap::ThreeMatrix<double> unit({1,0,0, 0,1,0, 0,0,1});

    SECTION( "Initialization" ) {
        yap::ThreeMatrix<double> m({});
        REQUIRE( m == zero);

        auto u = yap::unitMatrix<double, 3>();
        REQUIRE( u == unit);
    }

    SECTION( "Transpose" ) {
        yap::ThreeMatrix<double> m({1,2,3, 4,5,6, 7,8,9});
        yap::ThreeMatrix<double> m_T({1,4,7, 2,5,8, 3,6,9});

        REQUIRE( m == m_T );
    }

    SECTION( "+ -" ) {
        yap::ThreeMatrix<double> m({1,2,3, 4,5,6, 7,8,9});
        yap::ThreeMatrix<double> minus_m({-1,-2,-3, -4,-5,-6, -7,-8,-9});

        REQUIRE( m - m == zero);
        REQUIRE( -m == m_m );
        REQUIRE( -1. * m == m_m );
        REQUIRE( m + minus_m == zero );
        REQUIRE( m + m == 2 * m );
    }


}
