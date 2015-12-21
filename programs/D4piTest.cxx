#include "BreitWigner.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "FinalStateParticle.h"
#include "FourVector.h"
#include "InitialStateParticle.h"
#include "make_unique.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "SpinUtilities.h"
#include "WignerD.h"

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

#include <assert.h>
#include <iostream>
#include <memory>
#include <string>

//#include <callgrind.h>

#include "logging.h"

int main( int argc, char** argv)
{

    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    unsigned max2L(2 * 4);

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") : ".") + "/evt.pdl");

    // initial state particle
    double radialSize = 1.;
    std::shared_ptr<yap::InitialStateParticle> D = factory.createInitialStateParticle(421, radialSize);

    // final state particles
    std::shared_ptr<yap::FinalStateParticle> piPlus = factory.createFinalStateParticle(211);
    std::shared_ptr<yap::FinalStateParticle> piMinus = factory.createFinalStateParticle(-211);

    // Set final-state particles
    D->setFinalStateParticles({piPlus, piMinus, piPlus, piMinus});

    // rho rho
    std::shared_ptr<yap::Resonance> rho = factory.createResonance(113, radialSize, std::make_unique<yap::BreitWigner>());
    rho->addChannels(piPlus, piMinus, max2L);

    D->addChannels(rho, rho, max2L);

    // omega omega
    std::shared_ptr<yap::Resonance> omega = factory.createResonance(223, radialSize, std::make_unique<yap::BreitWigner>());
    omega->addChannels(piPlus, piMinus, max2L);

    D->addChannels(omega, omega, max2L);

    // rho omega
    D->addChannels(rho, omega, max2L);

    // a_1 channels
    std::shared_ptr<yap::Resonance> sigma = factory.createResonance(9000221, radialSize, std::make_unique<yap::BreitWigner>());
    sigma->addChannels(piPlus, piMinus, max2L);

    std::shared_ptr<yap::Resonance> a_1 = factory.createResonance(20213, radialSize, std::make_unique<yap::BreitWigner>());
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
    for (auto& pc : D->particleCombinations())
        std::cout << *pc << "\n";
    std::cout << "\n";

    std::cout << "\nFour momenta symmetrizations with " << D->fourMomenta().maxSymmetrizationIndex() + 1 << " indices \n";
    for (auto& pc : D->fourMomenta().particleCombinations())
        std::cout << *pc << ": " << D->fourMomenta().symmetrizationIndex(pc) << "\n";

    std::cout << "\nHelicity angles symmetrizations with " << D->helicityAngles().maxSymmetrizationIndex() + 1 << " indices \n";
    for (auto& pc : D->helicityAngles().particleCombinations())
        std::cout << *pc << ": " << D->helicityAngles().symmetrizationIndex(pc) << "\n";

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

    LOG(INFO) << "create dataPoints";

    for (unsigned int iEvt = 0; iEvt < 4; ++iEvt) {
        TGenPhaseSpace event;
        event.SetDecay(P, 4, masses);
        event.Generate();

        std::vector<yap::FourVector<double> > momenta;
        for (unsigned i = 0; i < 4; ++i) {
            TLorentzVector p = *event.GetDecay(i);
            momenta.push_back(yap::FourVector<double>({p.T(), p.X(), p.Y(), p.Z()}));

            DEBUG(yap::to_string(momenta.back()));
        }

        assert(D->addDataPoint(momenta));

        // test 4 momenta calculation
        DEBUG("calculated 4 momenta from invariant masses:");
        auto fourVecs = D->fourMomenta().calculateFourMomenta(D->dataSet().back());
        for (auto v : fourVecs)
            DEBUG(yap::to_string(v));
    }

    LOG(INFO) << "done creating dataPoints";

    D->dataSet()[0].printDataSize();

    LOG(INFO) << "Printing data:";
    for (unsigned d = 0; d < D->dataSet().size(); ++d) {
        LOG(INFO) << "  DataPoint " << d;
        for (auto& v : D->dataSet()[d].finalStateFourMomenta()) {
            LOG(INFO) << v[0] << "; " << v[1] << ", " << v[2] << ", " << v[3];
        }
    }


    // create data partitions
    D->setDataPartitions(yap::createDataPartitionsBlock(D->dataSet(), 2));

    // to test amplitude calculation, set all free amps to 1
    auto freeAmps = D->freeAmplitudes();

    LOG(INFO) << freeAmps.size() << " free amplitudes";

    for (auto& a : freeAmps)
        a->setValue(yap::Complex_1);

    //CALLGRIND_START_INSTRUMENTATION

    // do several loops over all dataPartitions
    for (unsigned i = 0; i < 2; ++i) {

        // change amplitudes
        for (auto& a : freeAmps)
            a->setValue(0.9 * a->value());
        DEBUG("===================================================================================================================== ");

        double logA(0);

        if (false) {
            // multi threaded
            logA = D->sumOfLogsOfSquaredAmplitudes();
        } else {
            // update global calculationStatuses before looping over partitions
            D->updateGlobalCalculationStatuses();

            // loop over partitions
            for (yap::DataPartitionBase* partition : D->dataPartitions()) {
                DEBUG("calculate logA for partition " << partition->index() << " ---------------------------------------------------------------------------");
                logA += D->partialSumOfLogsOfSquaredAmplitudes(partition);
            }

            // set parameter flags to unchanged after looping over all partitions
            D->setParameterFlagsToUnchanged();
        }

        LOG(INFO) << "logA = " << logA;
    }

    //CALLGRIND_STOP_INSTRUMENTATION

    /*
        for (auto& a : freeAmps)
            a->setValue(0.5 * a->value());

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);

        // only change some amps
        unsigned i(0);
        for (auto& a : freeAmps) {
            if (i++ % 2 == 0)
                a->setValue(1.);
        }

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);

        // set to zero
        for (auto& a : freeAmps) {
            a->setValue(0.);
        }

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);
    */


    std::cout << "alright! \n";
}
