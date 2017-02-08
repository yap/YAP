#include <catch.hpp>

#include <Attributes.h>
#include <BlattWeisskopf.h>
#include <DecayingParticle.h>
#include <HelicityFormalism.h>
#include <Model.h>

#include "helperFunctions.h"

TEST_CASE( "Attributes" )
{

    /// test throw on nullptr for attribute_of
    REQUIRE_THROWS_AS( yap::has_a_mass()((yap::DecayingParticle*)nullptr), yap::exceptions::Exception );

    auto M = d3pi<yap::HelicityFormalism>();
    auto rho = std::static_pointer_cast<yap::DecayingParticle>(particle(*M, yap::is_named("rho0")));
    auto piPlus = particle(*M, yap::is_named("pi+"));
    
    yap::orbital_angular_momentum get_l;
    REQUIRE(get_l(rho->blattWeisskopfs().at(1)) == 1);
    REQUIRE(get_l(rho->decayTrees()[0]) == 1);

    yap::spin_angular_momentum get_m;
    REQUIRE(get_m(rho) == 2);
    REQUIRE(get_m(rho->decayTrees()[0]) == 0);

    REQUIRE_THROWS_AS(yap::mass_parameter()(piPlus), yap::exceptions::Exception );
    
}
