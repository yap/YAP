#include <catch.hpp>

#include "helperFunctions.h"

#include <DataSet.h>
#include <Exceptions.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <Model.h>

TEST_CASE("DataSet")
{

    auto M3 = d3pi<yap::HelicityFormalism>();
    auto D3 = M3->createDataSet(1);

    // auto mass_axes = M->massAxes();
    // D.push_back(calculate_four_momenta(1.869, M->finalStateParticles(), mass_axes, std::vector<double>(mass_axes.size(), 1)));

    auto M4 = d4pi();
    auto D4 = M4->createDataSet(1);

    REQUIRE( D3.consistent(D4[0]) == false );

    REQUIRE_THROWS_AS( D3.push_back(D4[0]), yap::exceptions::Exception );
    REQUIRE_NOTHROW( D3.push_back(D3[0]) );

    auto d4 = D4[0];
    REQUIRE_THROWS_AS( D3.push_back(std::move(d4)), yap::exceptions::Exception );
    auto d3 = D3[0];
    REQUIRE_NOTHROW( D3.push_back(std::move(d3)) );

    REQUIRE_NOTHROW( D3.erase(D3.begin() + 1, D3.end()) );
    REQUIRE_THROWS_AS( D3.erase(D4.begin(), D3.end()), yap::exceptions::Exception );
}
