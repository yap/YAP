#include "BreitWigner.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "make_unique.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "SpinUtilities.h"
#include "WignerD.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <assert.h>
#include <iostream>
#include <string>

//#include <callgrind.h>

#include "logging.h"

int main( int argc, char** argv)
{

    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    unsigned max2L(2 * 4);

    yap::ParticleFactory factory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // initial state particle
    double radialSize = 1.;
    auto D = factory.createInitialStateParticle(421, radialSize);

    // final state particles
    auto piPlus = factory.createFinalStateParticle(211, {0, 2});
    auto piMinus = factory.createFinalStateParticle(-211, {1, 3});

    // rho rho
    auto rho = factory.createResonance(113, radialSize, std::make_unique<yap::BreitWigner>());
    rho->addChannels(piPlus, piMinus, max2L);

    D->addChannels(rho, rho, max2L);

    // omega omega
    auto omega = factory.createResonance(223, radialSize, std::make_unique<yap::BreitWigner>());
    omega->addChannels(piPlus, piMinus, max2L);

    D->addChannels(omega, omega, max2L);

    // rho omega
    D->addChannels(rho, omega, max2L);

    // a_1 channels
    auto sigma = factory.createResonance(9000221, radialSize, std::make_unique<yap::BreitWigner>());
    sigma->addChannels(piPlus, piMinus, max2L);

    auto a_1 = factory.createResonance(20213, radialSize, std::make_unique<yap::BreitWigner>());
    a_1->addChannels(sigma, piPlus, max2L);

    a_1->addChannels(rho, piPlus, max2L);

    D->addChannels(a_1, piMinus, max2L);


    // R pi pi channels
    //yap::Resonance* f_0_980 = factory.createResonanceBreitWigner(9000221, radialSize);
    //factory.createChannel(f_0_980, piPlus, piMinus, 0);


    // consistency and optimizations
    assert(D->prepare());
    std::cout << "consistent! \n";

    // print stuff
    //yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\n" << D->particleCombinations().size() << " D symmetrizations \n";
    /*for (auto& pc : D->particleCombinations())
        std::cout << std::string(*pc) << "\n";
    std::cout << "\n";*/

    std::cout << "\nFour momenta symmetrizations with " << D->fourMomenta().maxSymmetrizationIndex() + 1 << " indices \n";
    /*for (auto& pc : D->fourMomenta().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->fourMomenta().symmetrizationIndex(pc) << "\n";*/

    std::cout << "\nHelicity angles symmetrizations with " << D->helicityAngles().maxSymmetrizationIndex() + 1 << " indices \n";
    /*for (auto& pc : D->helicityAngles().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->helicityAngles().symmetrizationIndex(pc) << "\n";*/

    D->printDecayChain();
    std::cout << "\n";

    D->printSpinAmplitudes();
    D->printDataAccessors(false);
    //D->printDataAccessors();



    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    Double_t masses[4] = { piPlus->mass()->value(), piMinus->mass()->value(),
                           piPlus->mass()->value(), piMinus->mass()->value()
                         };

    LOG(INFO) << "create dataPoint";

    TGenPhaseSpace event;
    event.SetDecay(P, 4, masses);
    event.Generate();

    std::vector<TLorentzVector> momenta;
    for (unsigned i = 0; i < 4; ++i)
        momenta.push_back(*event.GetDecay(i));

    // Dalitz coordinates: m^2_{12}, m^2_{13}, m^2_{14}, m^2_{23}, m^2_{24}
    D->addDataPoint(momenta);

    yap::DataPoint& d = D->dataSet()[0];

    LOG(INFO) << "inv masses [GeV]";
    LOG(INFO) << std::string(*D->fourMomenta().InitialStatePC_) << ": " << D->fourMomenta().m(d, D->fourMomenta().InitialStatePC_);
    for (auto& pc : D->fourMomenta().FinalStatePC_)
        LOG(INFO) << std::string(*pc) << ": " << D->fourMomenta().m(d, pc);
    for (auto& pcV : D->fourMomenta().PairPC_)
        for (auto& pc : pcV)
            if (pc)
                LOG(INFO) << std::string(*pc) << ": " << D->fourMomenta().m(d, pc);

    LOG(INFO) << "inv mass squares [GeV^2]";
    LOG(INFO) << std::string(*D->fourMomenta().InitialStatePC_) << ": " << D->fourMomenta().m2(d, D->fourMomenta().InitialStatePC_);
    for (auto& pc : D->fourMomenta().FinalStatePC_)
        LOG(INFO) << std::string(*pc) << ": " << D->fourMomenta().m2(d, pc);
    for (auto& pcV : D->fourMomenta().PairPC_)
        for (auto& pc : pcV)
            if (pc)
                LOG(INFO) << std::string(*pc) << ": " << D->fourMomenta().m2(d, pc);

    // test formula
    double epsilon = D->fourMomenta().m2(d, D->fourMomenta().InitialStatePC_);
    for (auto& pc : D->fourMomenta().FinalStatePC_)
        epsilon -= (2.-4.) * D->fourMomenta().m2(d, pc);

    yap::ParticleCombinationVector pairPCs;
    unsigned i(0);
    for (auto& v : D->fourMomenta().PairPC_) {
        for (unsigned j=i++; j<v.size(); ++j)
            if (v[j])
                pairPCs.push_back(v[j]);
    }

    for (auto& pc : pairPCs)
        epsilon -= D->fourMomenta().m2(d, pc);

    LOG(INFO) << "epsilon = " << epsilon;

    // test
    unsigned symInd = D->fourMomenta().symmetrizationIndex(D->fourMomenta().PairPC_.at(0).at(1));
    LOG(INFO) << "symInd " << symInd;
    LOG(INFO) << "set value " << D->fourMomenta().M_.value(d, symInd);
    D->fourMomenta().M_.setValue(-1, d, symInd, 0u);

    for (auto& pc : D->fourMomenta().RecoilPC_)
        D->fourMomenta().M_.setValue(-1, d, D->fourMomenta().symmetrizationIndex(pc), 0u);


    // DO THE MAGIC
    D->fourMomenta().calculateMissingMasses(d);



    // result
    LOG(INFO) << "result:";

    LOG(INFO) << "inv mass squares [GeV^2]";
    LOG(INFO) << std::string(*D->fourMomenta().InitialStatePC_) << ": " << D->fourMomenta().m2(d, D->fourMomenta().InitialStatePC_);
    for (auto& pc : D->fourMomenta().FinalStatePC_)
        LOG(INFO) << std::string(*pc) << ": " << D->fourMomenta().m2(d, pc);
    for (auto& pcV : D->fourMomenta().PairPC_)
        for (auto& pc : pcV)
            if (pc)
                LOG(INFO) << std::string(*pc) << ": " << D->fourMomenta().m2(d, pc);

    epsilon = D->fourMomenta().m2(d, D->fourMomenta().InitialStatePC_);
    for (auto& pc : D->fourMomenta().FinalStatePC_)
        epsilon -= (2.-4.) * D->fourMomenta().m2(d, pc);
    for (auto& pc : pairPCs)
        epsilon -= D->fourMomenta().m2(d, pc);

    LOG(INFO) << "epsilon = " << epsilon;


    std::cout << "awlright! \n";
}
