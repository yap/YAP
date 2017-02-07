#include <catch.hpp>

#include <ParticleCombination.h>
#include <ParticleCombinationCache.h>

#include "helperFunctions.h"

#include <memory>

TEST_CASE( "ParticleCombination" )
{

    auto M = d4pi();

    for (const auto& pc : M->particleCombinationCache()) {
        if (!pc.lock())
            continue;
        REQUIRE_NOTHROW(yap::to_string(*pc.lock()));
        REQUIRE_NOTHROW(yap::to_string_with_parent(*pc.lock()));
    }
}
