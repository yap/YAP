#include <catch.hpp>

#include <VariableStatus.h>

TEST_CASE("VariableStatus")
{

    SECTION ( "String Output") {

        REQUIRE_NOTHROW( yap::to_string(yap::VariableStatus::changed) );
        REQUIRE_NOTHROW( yap::to_string(yap::VariableStatus::fixed) );
        REQUIRE_NOTHROW( yap::to_string(yap::VariableStatus::unchanged) );
        
    }
    
}
