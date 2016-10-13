// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_fit.h"
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
#include <ZemachFormalism.h>

#include <algorithm>
#include <random>

int main()
{
    yap::plainLogs(el::Level::Info);

    unsigned i_model = 0;

    // create model
    std::string model_name;
    bat_fit* m;
    double D_mass(0);

    switch (i_model) {
        case 0:
            model_name = "D3PI";
            m = new bat_fit(d3pi_fit(model_name + "_fit", yap_model<yap::ZemachFormalism>()));
            D_mass = 1.86961; // D+
            break;
        case 1:
            model_name = "DKKPI";
            m = new bat_fit(dkkpi_fit(model_name + "_fit", yap_model<yap::HelicityFormalism>()));
            D_mass = 1.86961; // D+
            break;
        case 2:
            model_name = "D4PI";
            m = new bat_fit(d4pi_fit("D4pi_fit"));
            D_mass = 1.8648400; // D0
            break;
        default:
            LOG(ERROR) << "No model loaded";
    }

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(D_mass, m->axes(), m->model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);

    std::generate_n(std::back_inserter(m->integralData()), 100000,
            std::bind(yap::phsp<std::mt19937>, std::cref(*m->model()), D_mass, m->axes(), m2r, g, std::numeric_limits<unsigned>::max()));
    m->integralPartitions() = yap::DataPartitionBlock::create(m->integralData(), 6);
    LOG(INFO) << "Created " << m->integralData().size() << " data points (" << (m->integralData().bytes() * 1.e-6) << " MB)";

    yap::ImportanceSampler::calculate(m->modelIntegral(), m->integralPartitions());

    double sum(0);
    for (const auto& b2_dtvi : m->modelIntegral().integrals()) {
        auto ff = fit_fractions(b2_dtvi.second);
        for (size_t i = 0; i < ff.size(); ++i) {
            LOG(INFO) << to_string(*b2_dtvi.second.decayTrees()[i]) << "\t" << ff[i].value()*100. << " %";
            sum += ff[i].value();
        }
    }
    LOG(INFO) << "Sum = " << sum*100 << " %";

    delete m;

    return 0;
}
