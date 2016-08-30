/// \file test_InvariantMassBinning.cxx
/// \brief Test the binning of the mass axis.
/// \author Paolo Di Giglio

#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>

#include <DataSet.h>
#include <DataPoint.h>
#include <InvariantMassBinning.h>
#include <FourMomenta.h>
#include <MassRange.h>
#include <Model.h>
#include <ParticleCombination.h>

#include "helperFunctions.h"
#include <iostream>
#include <memory>
#include <vector>

const std::vector<double> monospaced_partition(const yap::MassRange& range, std::size_t bins)
{
    const double bin_width = (range[1] - range[0]) / bins;
    // Reserve the space for the partition
    std::vector<double> p;
    p.reserve(bins + 1);

    for (std::size_t i = 0; i < bins + 1; ++i)
        p.push_back(i * bin_width);

    return p;
}

inline const std::vector<yap::MassRange> model_mass_range(std::shared_ptr<yap::Model> model)
{
    return yap::mass_range(
               model->massAxes(),
               model->initialStateParticles().begin()->first,
               model->finalStateParticles());
}

inline const bool mass_belongs_to_bin(double mass, int bin, const std::vector<double>& bin_low_edges)
{
    if (bin == static_cast<int>(bin_low_edges.size()) - 1)
        return bin_low_edges[bin] <= mass;
    if (bin == -1)
        return mass < bin_low_edges[0];

    return bin_low_edges[bin] <= mass && mass < bin_low_edges[bin + 1];
}

class InvariantMassBinningTest : public yap::InvariantMassBinning
{
public:
    InvariantMassBinningTest(yap::Model& m, const std::vector<double>& bins) :
        yap::InvariantMassBinning(m, bins)
    {}

    void addPC(std::shared_ptr<yap::ParticleCombination> pc)
    { yap::InvariantMassBinning::addParticleCombination(pc); }
};

TEST_CASE( "InvariantMassBinning" )
{
    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    // Create the model
    auto model = create_model();
    REQUIRE(model->consistent());

    // Get the allowed mass range for the model mass axes
    auto m_range = model_mass_range(model);
//    for (const auto& it : m_range)
//        std::cerr << it[0] << " " << it[1] << std::endl;

    // Number of bins
    const unsigned int nBins = 10;
    // Vector of the low edges of the partition
    auto bin_low_edges = monospaced_partition(*m_range.begin(), nBins);
    // Create the binning
    InvariantMassBinningTest binning(*model, bin_low_edges);

    // Loop over particle combinations -> indices
    for (const auto& pc_i : model->fourMomenta()->symmetrizationIndices()) {
        // Add particle combinations to the binning
        binning.addPC(pc_i.first);
    }

    // Generate the DataSet
    const unsigned int nPoints = 1000;
    auto data_set = generate_data(model, nPoints);

    // Calculate and check the bins
    for (const auto& pc_i : model->fourMomenta()->symmetrizationIndices()) {
        for (auto& dp : data_set) {
            // Calculate the bin
            binning.calculate(dp, data_set);
            // Check the bin
            auto mass = model->fourMomenta()->m(dp, pc_i.first);
            auto bin  = binning.bin(dp, pc_i.first);
            auto bin_is_correct = mass_belongs_to_bin(mass, bin, bin_low_edges);
            if (!bin_is_correct)
                std::cout << "bin = " << bin << " "
                          << bin_low_edges[bin] << " <= "
                          << mass << " < "
                          << bin_low_edges[bin + 1] << std::endl;

            REQUIRE(bin_is_correct);
        }
    }
}
