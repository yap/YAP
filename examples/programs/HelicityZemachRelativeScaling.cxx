#include <BreitWigner.h>
#include <container_utils.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayTree.h>
#include <Filters.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <Resonance.h>
#include <ZemachFormalism.h>

#include "../../test/helperFunctions.h"


std::shared_ptr<yap::Model> create_model(std::unique_ptr<yap::SpinAmplitudeCache> formalism, unsigned L)
{
    auto M = std::make_shared<yap::Model>(std::move(formalism));

    yap::ParticleFactory F = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    double radialSize = 1.;

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);
    auto D_mass = F["D+"].Mass;

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    // Set final-state particles
    M->setFinalState(piPlus, piMinus, piPlus);

    // mass range
    auto axes = M->massAxes({{0, 1}, {1, 2}});
    auto m2_pipi_range = yap::squared(yap::mass_range(D_mass, axes, M->finalStateParticles()));
    double center = sqrt(m2_pipi_range[0][1] + 0.5*m2_pipi_range[0][0]);
    double width = 0.025;//sqrt(m2_pipi_range[0][1] - m2_pipi_range[0][0]) *100.;

    std::cout << "range "<< m2_pipi_range[0][0] << " .. " << m2_pipi_range[0][1] << "\n";
    std::cout << "center "<< center << "; width = " << width << "\n";

    auto pipi = yap::Resonance::create("pipi", yap::QuantumNumbers(2*L, 0), 3., std::make_shared<yap::BreitWigner>(center, width));
    pipi->addChannel(piPlus, piMinus);

    // D's channels
    D->addChannel(pipi, piPlus);

    M->addInitialStateParticle(D);

    return M;
}

double integrate(std::unique_ptr<yap::SpinAmplitudeCache> formalism, unsigned L)
{
        const unsigned nPoints = 5000;

        auto M = create_model(std::move(formalism), L);
        auto data = generate_data(*M, nPoints);
        auto partitions = yap::DataPartitionBlock::create(data, 4);

        yap::ModelIntegral mi(*M);

        yap::ImportanceSampler::calculate(mi, partitions);
        double smartIntegral = integral(mi).value();
        //LOG(INFO) << "Integral = " << smartIntegral;

        return smartIntegral;
}

int main () {
    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);


    for (unsigned L = 0; L < 3; ++L) {
        double z = integrate(std::make_unique<yap::ZemachFormalism>(), L);
        double h = integrate(std::make_unique<yap::HelicityFormalism>(), L);

        LOG(INFO) << "L = " << L << "; Ratio z/h = " << z/h;
    }

    return 0;
}
