#include <catch.hpp>

#include "helperFunctions.h"

#include <ParticleTable.h>
#include <PDL.h>
#include <QuantumNumbers.h>

TEST_CASE( "ParticleTable" )
{

    auto T = yap::read_pdl_file(find_pdl_file());

    // check for nonexistent entry
    REQUIRE_THROWS_AS( T[1234567890], yap::exceptions::Exception );
    REQUIRE_THROWS_AS( T["qcdAxion"], yap::exceptions::Exception );
    
    // test empty name throw
    REQUIRE_THROWS_AS( yap::ParticleTableEntry(12345, "", yap::QuantumNumbers(0, 0), 0.), yap::exceptions::Exception);
    // test negative mass throw
    REQUIRE_THROWS_AS( yap::ParticleTableEntry(12345, "test", yap::QuantumNumbers(0, 0), -1), yap::exceptions::Exception);

    yap::ParticleTableEntry pte(311, "test", yap::QuantumNumbers(0, 0), 0.);

    // check for replacing entry
    REQUIRE_NOTHROW( T.insert(pte) );
    
    // check for nonexistent parameter throw
    REQUIRE_THROWS_AS( get_nth_element(pte, 1, "test"), yap::exceptions::Exception );

    yap::ParticleTable T2;

    REQUIRE_NOTHROW( T += T2 );
    
}
