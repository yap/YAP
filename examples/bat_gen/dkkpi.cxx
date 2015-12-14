// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "dkkpi.h"

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

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <complex>

// ---------------------------------------------------------
dkkpi::dkkpi(std::string name)
    : BCModel(name)
{
    yap::plainLogs(el::Level::Debug);

    yap::ParticleFactory factory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]
    // use only L up to 4
    unsigned max2L(2 * 4);

    // final state particles
    auto kPlus  = factory.createFinalStateParticle(+321);
    auto kMinus = factory.createFinalStateParticle(-321);
    auto piPlus = factory.createFinalStateParticle(+211);

    // initial state particle
    D_ = factory.createInitialStateParticle(factory.pdgCode("D+"), radialSize);
    D_->setFinalStateParticles({kPlus, kMinus, piPlus});

    // phi
    auto phi = std::make_shared<yap::Resonance>(factory.quantumNumbers("phi"), 1010.e-3, "phi", radialSize, std::make_unique<yap::BreitWigner>());
    //auto phi = std::make_shared<yap::Resonance>(factory.quantumNumbers("phi"), 1310.e-3, "phi", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(phi->massShape()).width()->setValue(40e-3);
    phi->addChannels(kPlus, kMinus, max2L);

    // X_2
    auto X_2 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_2"), 1.2, "X_2", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(X_2->massShape()).width()->setValue(80e-3);
    X_2->addChannels(piPlus, kMinus, max2L);

    // Add channels to D
    //D_->addChannels(phi,      piPlus, max2L);
    D_->addChannels(X_2,      kPlus,  max2L);
    // D_->addChannels(f_2,      piPlus, max2L);
    // D_->addChannels(f_0_980,  piPlus, max2L);
    // D_->addChannels(f_0_1370, piPlus, max2L);
    // D_->addChannels(f_0_1500, piPlus, max2L);
    // D_->addChannels(sigma,    piPlus, max2L);

    // D_->addChannels(f_0_500_100,  piPlus, max2L);
    // D_->addChannels(f_0_1500_100, piPlus, max2L);

    D_->prepare();


    std::cout << "\n" << D_->particleCombinations().size() << " D symmetrizations \n";
    for (auto& pc : D_->particleCombinations())
        std::cout << std::string(*pc) << "\n";
    std::cout << "\n";

    std::cout << "\nFour momenta symmetrizations with " << D_->fourMomenta().maxSymmetrizationIndex() + 1 << " indices \n";
    for (auto& pc : D_->fourMomenta().particleCombinations())
        std::cout << std::string(*pc) << ": " << D_->fourMomenta().symmetrizationIndex(pc) << "\n";

    std::cout << "\nHelicity angles symmetrizations with " << D_->helicityAngles().maxSymmetrizationIndex() + 1 << " indices \n";
    for (auto& pc : D_->helicityAngles().particleCombinations())
        std::cout << std::string(*pc) << ": " << D_->helicityAngles().symmetrizationIndex(pc) << "\n";

    D_->printDecayChain();
    std::cout << "\n";

    D_->printSpinAmplitudes();
    D_->printDataAccessors(false);


    std::vector<std::shared_ptr<yap::ComplexParameter> > freeAmps = D_->freeAmplitudes();
    for (unsigned i = 0; i < freeAmps.size(); ++i)
        freeAmps[i]->setValue(yap::Complex_1);

    // unsigned i = 0;
    // freeAmps[i++]->setValue(yap::Complex_1);
    // freeAmps[i++]->setValue(yap::Complex_i);

    // freeAmps[i++]->setValue(std::polar(1.,     0.)); // rho
    // freeAmps[i++]->setValue(std::polar(2.1, -123. * TMath::Pi() / 180.)); // f_2
    // freeAmps[i++]->setValue(std::polar(1.4,   12. * TMath::Pi() / 180.)); // f_0_980
    // freeAmps[i++]->setValue(std::polar(1.3,  -21. * TMath::Pi() / 180.)); // f_0_1370
    // freeAmps[i++]->setValue(std::polar(1.1,  -44. * TMath::Pi() / 180.)); // f_0_1500
    // freeAmps[i++]->setValue(std::polar(3.7,   -3. * TMath::Pi() / 180.)); // sigma

    bool b = D_->initializeForMonteCarloGeneration(MCMCGetNChains());
    std::cout << "success = " << b << std::endl;
    std::cout << "number of data partitions = " << D_->dataPartitions().size() << std::endl;

    // initialize with random point in phasespace
    TLorentzVector P(0., 0., 0., D_->mass()->value());
    std::vector<double> masses(D_->finalStateParticles().size(), -1);
    for (unsigned i = 0; i < masses.size(); ++i)
        masses[i] = D_->finalStateParticles()[i]->mass()->value();
    TGenPhaseSpace event;
    event.SetDecay(P, masses.size(), &masses[0]);
    std::vector<yap::FourVector<double> > momenta(masses.size());
    // Generate events
    for (unsigned i = 0; i < MCMCGetNChains(); ++i) {
        event.Generate();
        for (unsigned i = 0; i < masses.size(); ++i) {
            TLorentzVector p = *event.GetDecay(i);
            momenta[i] = {p.T(), p.X(), p.Y(), p.Z()};
        }
        D_->dataSet()[i].setFinalStateFourMomenta(momenta);
        D_->fourMomenta().calculate(D_->dataSet()[i]);
    }


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
dkkpi::~dkkpi()
{
    // destructor
}

// ---------------------------------------------------------
double dkkpi::LogLikelihood(const std::vector<double>& parameters)
{
    // if (!std::isfinite(LogAPrioriProbability(parameters)))
    //     return -std::numeric_limits<double>::infinity();

    unsigned c = MCMCGetCurrentChain();

    DEBUG("Set mass squares to " << parameters[0] << ", " << parameters[1]);

    if (! D_->fourMomenta().setSquaredMasses(D_->dataSet()[c], DalitzAxes_, parameters)) {
        return -std::numeric_limits<double>::infinity();
    }

    DEBUG("call calculate");
    D_->calculate(D_->dataSet()[c]);

    // D_->updateGlobalCalculationStatuses();

    double L =  D_->logOfSquaredAmplitude(D_->dataSet()[c], c);
    // double L = D_->partialSumOfLogsOfSquaredAmplitudes(D_->dataPartitions()[c]);

    // D_->setParameterFlagsToUnchanged();

    DEBUG("dkkpi::LogLikelihood = " << L);

    return L;
}

// ---------------------------------------------------------
double dkkpi::LogAPrioriProbability(const std::vector<double>& parameters)
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
