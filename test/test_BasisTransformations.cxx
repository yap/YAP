#include <catch.hpp>
#include <catch_capprox.hpp>

#include <AmplitudeBasis.h>
#include <Constants.h>
#include <logging.h>
#include <Matrix.h>

TEST_CASE( "BasisTransformations" )
{

    // transversity
    yap::basis::transversity<double> t1(std::polar(0.624, yap::rad(357.)),
                                        std::polar(0.157, yap::rad(120.)),
                                        std::polar(0.384, yap::rad(163.)));
    // dummy variances
    auto& t1_cov = t1.covariance();
    t1_cov[0][0][0][0] = 1.210;  t1_cov[0][0][0][1] = 0.101;
    t1_cov[0][0][1][0] = 0.218;  t1_cov[0][0][1][1] = 1.951;

    t1_cov[1][1][0][0] = 2.176;  t1_cov[1][1][0][1] = 0.153;
    t1_cov[1][1][1][0] = 0.156;  t1_cov[1][1][1][1] = 1.686;

    t1_cov[2][2][0][0] = 1.782;  t1_cov[2][2][0][1] = 0.278;
    t1_cov[2][2][1][0] = 0.279;  t1_cov[2][2][1][1] = 2.064;

    REQUIRE(t1.covariance()[0][0][0][0] != 0);

    auto c1 = yap::basis::canonical<double>(t1);
    auto t2 = yap::basis::transversity<double>(c1);

    REQUIRE( t1.longitudinal()  == Catch::Detail::CApprox(t2.longitudinal()) );
    REQUIRE( t1.parallel()      == Catch::Detail::CApprox(t2.parallel()) );
    REQUIRE( t1.perpendicular() == Catch::Detail::CApprox(t2.perpendicular()) );
    REQUIRE( t1.covariance()    == t2.covariance() );

    auto h1 = yap::basis::helicity<double>(t1);
    auto h2 = yap::basis::helicity<double>(c1);

    REQUIRE( h1.zero()       == Catch::Detail::CApprox(h2.zero() ) );
    REQUIRE( h1.plus()       == Catch::Detail::CApprox(h2.plus() ) );
    REQUIRE( h1.minus()      == Catch::Detail::CApprox(h2.minus()) );
    REQUIRE( h1.covariance() == h2.covariance() );

    auto t3 = yap::basis::transversity<double>(h1);

    REQUIRE( t3.longitudinal()  == Catch::Detail::CApprox(t1.longitudinal() ) );
    REQUIRE( t3.parallel()      == Catch::Detail::CApprox(t1.parallel()     ) );
    REQUIRE( t3.perpendicular() == Catch::Detail::CApprox(t1.perpendicular()) );
    REQUIRE( t3.covariance()    == t1.covariance() );

    auto c2 = yap::basis::canonical<double>(h1);

    REQUIRE( c2[0]           == Catch::Detail::CApprox(c1[0]) );
    REQUIRE( c2[1]           == Catch::Detail::CApprox(c1[1]) );
    REQUIRE( c2[2]           == Catch::Detail::CApprox(c1[2]) );
    REQUIRE( c2.covariance() == c1.covariance());
}
