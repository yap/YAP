#include <catch.hpp>

#include <DataSet.h>
#include <Exceptions.h>
#include <HelicityFormalism.h>
#include <make_unique.h>
#include <Model.h>
#include <SpinAmplitudeCache.h>

#include "helperFunctions.h"

TEST_CASE( "Model" )
{

    SECTION( "SpinAmplitudeCache" ) {

        // test nullptr SAC
        REQUIRE_THROWS_AS( yap::Model(std::unique_ptr<yap::SpinAmplitudeCache>(nullptr)), yap::exceptions::Exception );

        // // test borrowed SAC
        // auto M = d4pi();
        // REQUIRE_THROWS( yap::Model(std::unique_ptr<yap::SpinAmplitudeCache>(M->spinAmplitudeCache())), yap::exceptions::Exception );
        
    }

    SECTION( "Intensity Throws" ) {

        yap::Model M(std::make_unique<yap::HelicityFormalism>());
        M.lock();
        yap::DataSet D1(M);
        yap::DataSet D2(M);
        
        REQUIRE_THROWS_AS( sum_of_log_intensity(M, D1), yap::exceptions::Exception );

        yap::DataPartitionVector DP;
        REQUIRE_THROWS_AS( sum_of_log_intensity(M, DP), yap::exceptions::Exception );

        DP.push_back(nullptr);
        REQUIRE_THROWS_AS( sum_of_log_intensity(M, DP), yap::exceptions::Exception );

        yap::DataPartitionVector DP2 = {&D1, &D2};
        REQUIRE_THROWS_AS( sum_of_log_intensity(M, DP2), yap::exceptions::Exception );

        // Need to implement in code!
        // yap::DataPartitionVector DP3 = {&D1, &D2};
        // REQUIRE_THROWS_AS( sum_of_log_intensity(M, DP3), yap::exceptions::Exception );

    }

    SECTION( "FSP / ISP" ) {

        auto M = d4pi();
        REQUIRE_THROWS_AS( M->setFinalState(M->finalStateParticles()), yap::exceptions::Exception );

        yap::Model N(std::make_unique<yap::HelicityFormalism>());
        REQUIRE_THROWS_AS( N.setFinalState(M->finalStateParticles()), yap::exceptions::Exception );

        REQUIRE_THROWS_AS( M->addInitialState(std::shared_ptr<yap::DecayingParticle>(nullptr)), yap::exceptions::Exception );

        M->lock();
        REQUIRE_THROWS_AS( M->addInitialState(M->initialStates()[0]), yap::exceptions::Exception );
        REQUIRE_THROWS_AS( N.addInitialState(M->initialStates()[0]), yap::exceptions::Exception );
        
    }
    
}
