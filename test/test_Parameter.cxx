#include <catch.hpp>

#include <Exceptions.h>
#include <Parameter.h>
#include <VariableStatus.h>

#include <complex>
#include <memory>
#include <vector>

TEST_CASE( "Parameter" )
{

    SECTION( "fixing" ) {

        yap::RealParameter P(0);

        // check that fixed parameter throws if changed
        P.variableStatus() = yap::VariableStatus::fixed;
        REQUIRE_THROWS_AS( P = 1, yap::exceptions::ParameterIsFixed );
        
    }

    SECTION( "RealParameter" ) {

        yap::RealParameter P(0);

        // check size
        REQUIRE( P.size() == 1 );

        // check setting by vector
        std::vector<double> X{1, 2, 3};
        P.setValue(X);
        REQUIRE( P.value() == X[0] );
        
    }

    SECTION( "NonnegativeRealParameter" ) {

        yap::NonnegativeRealParameter P(1);

        // check setting negative number
        REQUIRE_THROWS_AS( P = -1, yap::exceptions::Exception );
        REQUIRE_THROWS_AS( yap::NonnegativeRealParameter(-1), yap::exceptions::Exception );
        
    }

    SECTION( "PositiveRealParameter" ) {

        yap::PositiveRealParameter P(1);

        // check setting negative number
        REQUIRE_THROWS_AS( P = 0, yap::exceptions::Exception );
        REQUIRE_THROWS_AS( yap::PositiveRealParameter(0), yap::exceptions::Exception );
        
    }

    SECTION( "ComplexParameter" ) {

        yap::ComplexParameter P(0);

        // check size
        REQUIRE( P.size() == 2 );

        // check setting by vector
        std::vector<double> X{1, 2, 3};
        P.setValue(X);
        REQUIRE( P.value() == std::complex<double>(X[0], X[1]) );
        
    }

    SECTION( "ComplexComponentParameter" ) {

        auto P = std::make_shared<yap::ComplexParameter>(std::polar(12., 1.));
        yap::RealComponentParameter R(P);
        yap::ImaginaryComponentParameter I(P);

        REQUIRE( R.value() == real(P->value()) );
        REQUIRE( I.value() == imag(P->value()) );

        R = 5;
        I = 6;
        REQUIRE( P->value() == std::complex<double>(5, 6) );

        *P = std::complex<double>(15, 16);
        REQUIRE( R.value() == 15 );
        REQUIRE( I.value() == 16 );
    }
    
}
