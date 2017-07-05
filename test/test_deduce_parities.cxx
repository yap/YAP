#include <catch.hpp>

#include "helperFunctions.h"

#include <ParticleTable.h>
#include <PDL.h>

#include <deduce_parities.h>

TEST_CASE( "deduce_parities" )
{
    auto T = yap::read_pdl_file(find_pdl_file());

    REQUIRE_NOTHROW( deduce_meson_parities(T) );

    REQUIRE( T["pi+"].quantumNumbers().P() == -1 );
    REQUIRE( T["pi0"].quantumNumbers().P() == -1 );
    REQUIRE( T["pi-"].quantumNumbers().P() == -1 );

    REQUIRE( T["eta"].quantumNumbers().P() == -1 );

    REQUIRE( T["rho+"].quantumNumbers().P() == -1 );
    REQUIRE( T["rho0"].quantumNumbers().P() == -1 );
    REQUIRE( T["rho-"].quantumNumbers().P() == -1 );

    REQUIRE( T["K_S0"].quantumNumbers().P() == -1 );
    REQUIRE( T["K_L0"].quantumNumbers().P() == -1 );

    REQUIRE( T["omega"].quantumNumbers().P() == -1 );

    REQUIRE( T["phi"].quantumNumbers().P() == -1 );
    
    REQUIRE( T["a_0+"].quantumNumbers().P() == +1 );
    REQUIRE( T["a_00"].quantumNumbers().P() == +1 );
    REQUIRE( T["a_0-"].quantumNumbers().P() == +1 );

    REQUIRE( T["D+"].quantumNumbers().P() == -1 );
    REQUIRE( T["D0"].quantumNumbers().P() == -1 );
    REQUIRE( T["anti-D0"].quantumNumbers().P() == -1 );
    REQUIRE( T["D-"].quantumNumbers().P() == -1 );

    REQUIRE( T["B+"].quantumNumbers().P() == -1 );
    REQUIRE( T["B0"].quantumNumbers().P() == -1 );
    REQUIRE( T["anti-B0"].quantumNumbers().P() == -1 );
    REQUIRE( T["B-"].quantumNumbers().P() == -1 );

    // should fail to deduce f_0 parity
    REQUIRE( T["f_0"].quantumNumbers().P() == 0 );
}
