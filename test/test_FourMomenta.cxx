#include <catch.hpp>

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

    yap::Model M(std::make_unique<yap::HelicityFormalism>());
    auto D = M.createDataSet();

    REQUIRE_THROWS_AS( D.createDataPoint(std::vector<yap::FourVector<double> >()), yap::exceptions::EmptyFourMomentaVector );
    REQUIRE_THROWS_AS( D.createDataPoint({yap::FourVector<double>({1, 1, 1, 1})}), yap::exceptions::Exception );

}
