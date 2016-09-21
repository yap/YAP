// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_fit.h"
//#include "hist.h"
#include "models/d3pi.h"
#include "models/d4pi.h"
#include "models/dkkpi.h"
#include "tools.h"

#include <HelicityFormalism.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <PHSP.h>
#include <RelativisticBreitWigner.h>
#include <ZemachFormalism.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

int main()
{
    yap::plainLogs(el::Level::Info);


    // create model
    // auto m = dkkpi_fit(model_name + "_fit", yap_model<yap::HelicityFormalism>(), find_mass_axes(*t_pars));
    // auto m = d3pi_fit(model_name + "_fit", yap_model<yap::ZemachFormalism>(), find_mass_axes(*t_pars));
    auto m = d4pi_fit("D4pi_fit");

    //double D_mass = 1.86961;
    double D_mass = 1.8648400;

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(D_mass, m.axes(), m.model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    if (false) {
        m.integrationPointGenerator() = std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), D_mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max());
        // m.setNIntegrationPoints(4e4, 4e4);
        m.setNIntegrationPoints(20e3, 20e3, 6);
        LOG(INFO) << "Generating integration points on the fly";
    } else {
        // get FSP mass ranges
        auto m2r = yap::squared(mass_range(D_mass, m.axes(), m.model()->finalStateParticles()));
        
        // generate integration data
        std::generate_n(std::back_inserter(m.integralData()), 10000,
                        std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), D_mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max()));
        m.integralPartitions() = yap::DataPartitionBlock::create(m.integralData(), 6);
        LOG(INFO) << "Created " << m.integralData().size() << " data points (" << (m.integralData().bytes() * 1.e-6) << " MB)";

        yap::ImportanceSampler::calculate(m.modelIntegral(), m.integralPartitions());

        for (const auto& b2_dtvi : m.modelIntegral().integrals()) {
            auto ff = fit_fractions(b2_dtvi.second);
            for (size_t i = 0; i < ff.size(); ++i) {
                LOG(INFO) << to_string(*b2_dtvi.second.decayTrees()[i]) << "\t" << ff[i].value()*100. << " %";            }
        }
    }


    return 0;
}
