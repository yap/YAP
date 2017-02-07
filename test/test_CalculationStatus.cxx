#include <catch.hpp>

#include <CalculationStatus.h>

TEST_CASE("CalculationStatus")
{

    SECTION ( "String Output") {

        REQUIRE_NOTHROW( yap::to_string(yap::CalculationStatus::calculated) );
        REQUIRE_NOTHROW( yap::to_string(yap::CalculationStatus::uncalculated) );
        
    }
    
}
