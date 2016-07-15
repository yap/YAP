#include <catch.hpp>
#include <catch_capprox.hpp>

#include <AmplitudeBasis.h>
#include <Constants.h>
#include <logging.h>

TEST_CASE( "BasisTransformations" )
{

    // transversity
    yap::basis::transversity<double> t1(std::polar(0.624, yap::rad(357.)),
                                        std::polar(0.157, yap::rad(120.)),
                                        std::polar(0.384, yap::rad(163.)));

    auto c1 = yap::basis::canonical<double>(t1);
    auto t2 = yap::basis::transversity<double>(c1);

    REQUIRE( t1.longitudinal  == Catch::Detail::CApprox(t2.longitudinal) );
    REQUIRE( t1.parallel      == Catch::Detail::CApprox(t2.parallel) );
    REQUIRE( t1.perpendicular == Catch::Detail::CApprox(t2.perpendicular) );

    auto h1 = yap::basis::helicity<double>(t1);
    auto h2 = yap::basis::helicity<double>(c1);

    REQUIRE( h1.zero  == Catch::Detail::CApprox(h2.zero ) );
    REQUIRE( h1.plus  == Catch::Detail::CApprox(h2.plus ) );
    REQUIRE( h1.minus == Catch::Detail::CApprox(h2.minus) );

    auto t3 = yap::basis::transversity<double>(h1);

    REQUIRE( t3.longitudinal  == Catch::Detail::CApprox(t1.longitudinal ) );
    REQUIRE( t3.parallel      == Catch::Detail::CApprox(t1.parallel     ) );
    REQUIRE( t3.perpendicular == Catch::Detail::CApprox(t1.perpendicular) );

    auto c2 = yap::basis::canonical<double>(h1);

    REQUIRE( c2[0] == Catch::Detail::CApprox(c1[0]) );
    REQUIRE( c2[1] == Catch::Detail::CApprox(c1[1]) );
    REQUIRE( c2[2] == Catch::Detail::CApprox(c1[2]) );
}
