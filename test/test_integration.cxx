#include <catch.hpp>

#include "helperFunctions.h"

#include <CompensatedSum.h>
#include <DataPartition.h>
#include <DataPoint.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <Model.h>
#include <ModelIntegral.h>

#include <future>
#include <memory>
#include <vector>

/**
 *  Test the integration
 */

TEST_CASE("integration")
{
    // disable debug logs in test
    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    const unsigned nPoints = 1000;

    auto M = d4pi();
    auto data = generate_data(*M, nPoints);
    auto partitions = yap::DataPartitionBlock::create(data, 4);

    // integrating lambda
    auto sum_of_intensities = [M](yap::DataPartition& D)
        {
            M->calculate(D);
            yap::CompensatedSum<double> I(0.);
            for (auto& d : D)
                I += intensity(*M, d);
            return I.sum;
        };
    
    // create threads for calculation on each partition
    std::vector<std::future<double> > partial_sums;
    partial_sums.reserve(partitions.size());
    for (auto& P : partitions)
        partial_sums.push_back(std::async(std::launch::async, sum_of_intensities, std::ref(*P)));

    // wait for each partition to finish calculating
    double bruteForceIntegral = std::accumulate(partial_sums.begin(), partial_sums.end(), 0., [](double l, auto& s){return l + s.get();});

    DEBUG("bruteForceIntegral = " << bruteForceIntegral);

    yap::ModelIntegral mi(*M);

    yap::ImportanceSampler::calculate(mi, partitions);
    double smartIntegral = integral(mi).value();
    DEBUG("smartIntegral = " << smartIntegral);

    REQUIRE(bruteForceIntegral / nPoints == Approx(smartIntegral));
}

