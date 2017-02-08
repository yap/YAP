#include <catch.hpp>

#include <Attributes.h>
#include <DecayingParticle.h>

TEST_CASE( "Attributes" )
{

    /// test throw on nullptr for attribute_of
    REQUIRE_THROWS_AS( yap::has_a_mass()((yap::DecayingParticle*)nullptr), yap::exceptions::Exception );

}
