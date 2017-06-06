#include <catch.hpp>

#include "helperFunctions.h"

#include <DataSet.h>
#include <Exceptions.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <HelicityFormalism.h>
#include <make_unique.h>
#include <Model.h>

#include <vector>

TEST_CASE( "FourMomenta" )
{

    auto M = d3pi<yap::HelicityFormalism>();

    REQUIRE(M->fourMomenta()->consistent());
    
    auto D = M->createDataSet();

    REQUIRE_THROWS_AS( D.createDataPoint(std::vector<yap::FourVector<double> >()), yap::exceptions::EmptyFourMomentaVector );
    REQUIRE_THROWS_AS( D.createDataPoint({yap::FourVector<double>({1, 1, 1, 1})}), yap::exceptions::Exception );

    auto mass_axes = M->massAxes();
    // empty FSPs
    REQUIRE_THROWS_AS( calculate_four_momenta(1.869, yap::FinalStateParticleVector(), mass_axes, std::vector<double>(mass_axes.size(), 1)),
                       yap::exceptions::Exception );
    // too many m2
    REQUIRE_THROWS_AS( calculate_four_momenta(1.869, M->finalStateParticles(), mass_axes, std::vector<double>(mass_axes.size() + 1, 1)),
                       yap::exceptions::Exception );
    // negative initial mass
    REQUIRE_THROWS_AS( calculate_four_momenta(-1.869, M->finalStateParticles(), mass_axes, std::vector<double>(mass_axes.size(), 1)),
                       yap::exceptions::Exception );
    // negative mass
    REQUIRE_THROWS_AS( calculate_four_momenta(1.869, M->finalStateParticles(), mass_axes, std::vector<double>(mass_axes.size(), -1)),
                       yap::exceptions::Exception );
    
    D.push_back(calculate_four_momenta(1.869, M->finalStateParticles(), mass_axes, std::vector<double>(mass_axes.size(), 1)));
    REQUIRE_NOTHROW( M->fourMomenta()->massesString(D[0]) );
}
