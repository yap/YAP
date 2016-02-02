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

    // final state particles
    std::shared_ptr<yap::FinalStateParticle> kPlus  = factory.createFinalStateParticle(+321);
    std::shared_ptr<yap::FinalStateParticle> kMinus = factory.createFinalStateParticle(-321);
    std::shared_ptr<yap::FinalStateParticle> piPlus = factory.createFinalStateParticle(+211);

    // initial state particle
    D_ = factory.createInitialStateParticle(factory.pdgCode("D+"), radialSize);
    D_->setFinalStateParticles({kPlus, kMinus, piPlus});

    // phi
    std::shared_ptr<yap::Resonance> phi = std::make_shared<yap::Resonance>(factory.quantumNumbers("phi"), 1010.e-3, "phi", radialSize, std::make_unique<yap::BreitWigner>());
    //std::shared_ptr<yap::Resonance> phi = std::make_shared<yap::Resonance>(factory.quantumNumbers("phi"), 1310.e-3, "phi", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(phi->massShape()).width()->setValue(40e-3);
    phi->addChannel({kPlus, kMinus});

    // X_2
    std::shared_ptr<yap::Resonance> X_2 = std::make_shared<yap::Resonance>(factory.quantumNumbers("f_2"), 1.2, "X_2", radialSize, std::make_unique<yap::BreitWigner>());
    static_cast<yap::BreitWigner&>(X_2->massShape()).width()->setValue(80e-3);
    X_2->addChannel({piPlus, kMinus});

    // Add channels to D
    D_->addChannel({phi, piPlus});
    D_->addChannel({X_2, kPlus});

    D_->prepare();

    std::vector<std::shared_ptr<yap::ComplexParameter> > freeAmps = D_->freeAmplitudes();
    for (unsigned i = 0; i < freeAmps.size(); ++i)
        freeAmps[i]->setValue(yap::Complex_1);

    D_->initializeForMonteCarloGeneration(3);
    SetNChains(3);
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
    for (unsigned i = 0; i < GetNChains(); ++i) {
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

    unsigned c = GetCurrentChain();

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
