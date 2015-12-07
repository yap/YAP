// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "d3pi.h"

#include <BAT/BCMath.h>

#include "BreitWigner.h"
#include "Constants.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "make_unique.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"

#include <complex>

// ---------------------------------------------------------
d3pi::d3pi(std::string name)
    : BCModel(name)
{
    yap::ParticleFactory factory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]
    // use only L up to 4
    unsigned max2L(2 * 4);

    // final state particles
    auto piPlus = factory.createFinalStateParticle(211);
    auto piMinus = factory.createFinalStateParticle(-211);

    // initial state particle
    D_ = factory.createInitialStateParticle(factory.pdgCode("D+"), radialSize);
    D_->setFinalStateParticles({piPlus, piMinus, piPlus});

    // rho
    auto rho = std::make_shared<yap::Resonance>(factory.quantumNumbers("rho0"), 0.775, "rho", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(rho->massShape()).width()->setValue(0.149);
    rho->addChannels(piPlus, piMinus, max2L);

    // // f_2(1270)
    auto f_2 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_2"), 1.275, "f_2", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_2->massShape()).width()->setValue(0.185);
    f_2->addChannels(piPlus, piMinus, max2L);

    // f_0(980)
    auto f_0_980 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 0.980, "f_0_980", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_980->massShape()).width()->setValue(0.329);
    f_0_980->addChannels(piPlus, piMinus, max2L);

    // f_0(1370)
    auto f_0_1370 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_1370->massShape()).width()->setValue(0.250);
    f_0_1370->addChannels(piPlus, piMinus, max2L);

    // f_0(1500)
    auto f_0_1500 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.507, "f_0_1500", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_1500->massShape()).width()->setValue(0.109);
    f_0_1500->addChannels(piPlus, piMinus, max2L);

    // sigma a.k.a. f_0(500)
    auto sigma = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 0.800, "sigma", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(sigma->massShape()).width()->setValue(0.800);
    sigma->addChannels(piPlus, piMinus, max2L);

    // f_0(500+100i)
    auto f_0_500_100 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), .500, "f_0_500_100", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_500_100->massShape()).width()->setValue(0.100);
    f_0_500_100->addChannels(piPlus, piMinus, max2L);

    // f_0(1500+100i)
    auto f_0_1500_100 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_0"), 1.500, "f_0_1500_100", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(f_0_1500_100->massShape()).width()->setValue(0.100);
    f_0_1500_100->addChannels(piPlus, piMinus, max2L);


    // Add channels to D
    D_->addChannels(rho,      piPlus, max2L);
    D_->addChannels(f_2,      piPlus, max2L);
    // D_->addChannels(f_0_980,  piPlus, max2L);
    // D_->addChannels(f_0_1370, piPlus, max2L);
    // D_->addChannels(f_0_1500, piPlus, max2L);
    // D_->addChannels(sigma,    piPlus, max2L);

    // D_->addChannels(f_0_500_100,  piPlus, max2L);
    // D_->addChannels(f_0_1500_100, piPlus, max2L);

    D_->prepare();

    std::vector<std::shared_ptr<yap::ComplexParameter> > freeAmps = D_->freeAmplitudes();
    for (unsigned i = 0; i < freeAmps.size(); ++i)
        freeAmps[i]->setValue(yap::Complex_1);

    // unsigned i = 0;
    // freeAmps[i++]->setValue(std::polar(1.,     0.)); // rho
    // freeAmps[i++]->setValue(std::polar(2.1, -123. * TMath::Pi() / 180.)); // f_2
    // freeAmps[i++]->setValue(std::polar(1.4,   12. * TMath::Pi() / 180.)); // f_0_980
    // freeAmps[i++]->setValue(std::polar(1.3,  -21. * TMath::Pi() / 180.)); // f_0_1370
    // freeAmps[i++]->setValue(std::polar(1.1,  -44. * TMath::Pi() / 180.)); // f_0_1500
    // freeAmps[i++]->setValue(std::polar(3.7,   -3. * TMath::Pi() / 180.)); // sigma

    bool b = D_->initializeForMonteCarloGeneration(3);
    std::cout << "success = " << b << std::endl;
    std::cout << "number of data partitions = " << D_->dataPartitions().size() << std::endl;

    DalitzAxes_ = D_->fourMomenta().getDalitzAxes({{0, 1}, {1, 2}});

    for (auto& pc : DalitzAxes_) {
        std::string axis_label = "m2_";
        for (auto& d : pc->daughters())
            axis_label += std::to_string(d->indices()[0]);
        auto mrange = D_->getMassRange(pc);
        AddParameter(axis_label, pow(mrange[0], 2), pow(mrange[1], 2));
        std::cout << "Added parameter " << axis_label
                  << " with range = [" << pow(mrange[0], 2) << ", " << pow(mrange[1], 2) << "]"
                  << std::endl;
    }

    m2_P = pow(D_->mass()->value(), 2);
    m2_a = pow(D_->finalStateParticles()[0]->mass()->value(), 2);
    m2_b = pow(D_->finalStateParticles()[1]->mass()->value(), 2);
    m2_c = pow(D_->finalStateParticles()[2]->mass()->value(), 2);

    // Define parameters here in the constructor. For example:
    // AddParameter("mu",-2,1,"#mu");
    // And set priors, if using built-in priors. For example:
    // GetParamater("mu").SetPrior(new BCPriorGaus(-1, 0.25));
}

// ---------------------------------------------------------
d3pi::~d3pi()
{
    // destructor
}

// ---------------------------------------------------------
double d3pi::LogLikelihood(const std::vector<double>& parameters)
{
    // if (!std::isfinite(LogAPrioriProbability(parameters)))
    //     return -std::numeric_limits<double>::infinity();

    unsigned c = MCMCGetCurrentChain();

    D_->fourMomenta().setSquaredMasses(D_->dataSet()[c], DalitzAxes_, parameters);

    // D_->updateGlobalCalculationStatuses();

    double L =  D_->logOfSquaredAmplitude(D_->dataSet()[c], c);
    // double L = D_->partialSumOfLogsOfSquaredAmplitudes(D_->dataPartitions()[c]);

    // D_->setParameterFlagsToUnchanged();

    return L;
}

// ---------------------------------------------------------
double d3pi::LogAPrioriProbability(const std::vector<double>& parameters)
{
    double m2_ab = parameters[0];
    double m2_bc = parameters[1];

    if (m2_ab < 0 or m2_bc < 0)
        return 0;

    if (m2_b < 0)
        m2_b = m2_a;
    if (m2_c < 0)
        m2_c = m2_a;

    if (m2_ab < m2_a + m2_b + 2 * sqrt(m2_a * m2_b) or m2_ab > m2_P + m2_c - 2 * sqrt(m2_P * m2_c))
        return 0;

    double Eb = (m2_ab - m2_a + m2_b) / 2 / sqrt(m2_ab);
    double Ec = (m2_P - m2_ab - m2_c) / 2 / sqrt(m2_ab);
    double Pb = sqrt(Eb * Eb - m2_b);
    double Pc = sqrt(Ec * Ec - m2_c);

    if (fabs(m2_bc - m2_b - m2_c - 2 * Eb * Ec) <= 2 * Pb * Pc)
        return 0;

    return -std::numeric_limits<double>::infinity();
}
