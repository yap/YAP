#include "BreitWigner.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "DataSet.h"
#include "DecayTree.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "FourVector.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "make_unique.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "SpinAmplitudeCache.h"
#include "Resonance.h"
#include "WignerD.h"

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <assert.h>
#include <iostream>
#include <memory>
#include <string>

//#include <callgrind.h>

int main( int argc, char** argv)
{

    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    yap::Model M(std::make_unique<yap::HelicityFormalism>());

//    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    yap::ParticleFactory factory("../../../data/evt.pdl");

    double radialSize = 1.;

    // initial state particle
    auto D = factory.decayingParticle(421, radialSize);

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // Set final-state particles
    M.setFinalState({piPlus, piMinus, piPlus, piMinus});

    // sigma
    auto sigma = factory.resonance(9000221, radialSize, std::make_shared<yap::BreitWigner>());
    sigma->addChannel({piPlus, piMinus});

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::BreitWigner>());
    rho->addChannel({piPlus, piMinus});

    // omega
    auto omega = factory.resonance(223, radialSize, std::make_shared<yap::BreitWigner>());
    omega->addChannel({piPlus, piMinus});

    // a_1
    auto a_1 = factory.resonance(20213, radialSize, std::make_shared<yap::BreitWigner>());
    a_1->addChannel({sigma, piPlus});
    a_1->addChannel({rho,   piPlus});

    // D's channels
    D->addChannel({rho, rho});
    D->addChannel({omega, omega});
    D->addChannel({rho, omega});
    D->addChannel({a_1, piMinus});

    // R pi pi channels
    //yap::Resonance* f_0_980 = factory.resonanceBreitWigner(9000221, radialSize);
    //factory.createChannel(f_0_980, piPlus, piMinus, 0);

    // check consistency
    if (M.consistent())
        LOG(INFO) << "consistent!";
    else
        LOG(INFO) << "inconsistent!";

    // print stuff
    //yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\n" << D->particleCombinations().size() << " D symmetrizations \n";
    for (auto& pc : D->particleCombinations())
        std::cout << *pc << "\n";
    std::cout << "\n";

    std::cout << "\nFour momenta symmetrizations with " << M.fourMomenta()->nSymmetrizationIndices() << " indices \n";
    for (auto& pc : M.fourMomenta()->particleCombinations())
        std::cout << *pc << ": " << M.fourMomenta()->symmetrizationIndex(pc) << "\n";

    std::cout << "\nHelicity angles symmetrizations with " << M.helicityAngles()->nSymmetrizationIndices() << " indices \n";
    for (auto& pc : M.helicityAngles()->particleCombinations())
        std::cout << *pc << ": " << M.helicityAngles()->symmetrizationIndex(pc) << "\n";

    D->printDecayChain();
    std::cout << "\n";

    std::cout << *M.spinAmplitudeCache() << std::endl;
    M.printDataAccessors(false);

    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    Double_t masses[4] = { piPlus->mass()->value(), piMinus->mass()->value(), piPlus->mass()->value(), piMinus->mass()->value() };

    LOG(INFO) << "create dataPoints";

    // create data set
    unsigned nPoints = 1;
    yap::DataSet data = M.createDataSet(nPoints);
    yap::DataSet dataTest = M.createDataSet(nPoints);

    for (unsigned int iEvt = 0; iEvt < nPoints; ++iEvt) {
        TGenPhaseSpace event;
        event.SetDecay(P, 4, masses);
        event.Generate();

        std::vector<yap::FourVector<double> > momenta;
        for (unsigned i = 0; i < 4; ++i) {
            TLorentzVector p = *event.GetDecay(i);
            momenta.push_back(yap::FourVector<double>({p.T(), p.X(), p.Y(), p.Z()}));

            DEBUG(yap::to_string(momenta.back()));
        }

        data[iEvt].setFinalStateMomenta(momenta);
        dataTest[iEvt].setFinalStateMomenta(momenta);

    }

    LOG(INFO) << "done creating dataPoints";

    LOG(INFO) << data[0].dataSizeString();

    LOG(INFO) << "Printing data:";
    for (unsigned d = 0; d < data.points().size(); ++d) {
        LOG(INFO) << "  DataPoint " << d;
        for (auto& v : M.fourMomenta()->finalStateMomenta(data[d]))
            LOG(INFO) << yap::to_string(v);
    }

    // create data partitions
    unsigned nChains = 1;
    auto parts = yap::DataPartitionWeave::create(data, nChains);
    //auto partsTest = yap::DataPartitionWeave::create(dataTest, nChains);

    for (auto& dt : M.decayTrees())
        LOG(INFO) << to_string(dt);

    auto freeAmps = M.freeAmplitudes();

    LOG(INFO) << freeAmps.size() << " free amplitudes";

    //CALLGRIND_START_INSTRUMENTATION

    // do several loops over all dataPartitions
    for (unsigned i = 0; i < 100; ++i) {

        // change amplitudes
        if (gRandom->Uniform() > 0.5) {
            for (auto& a : freeAmps) {
                if (a->variableStatus() != yap::VariableStatus::fixed and gRandom->Uniform() > 0.5)
                    a->setValue(gRandom->Uniform(0.95, 1.052631579) * a->value());
            }
        }
        DEBUG("===================================================================================================================== ");

        std::cout << "Variable status after changing:    \n";
        M.printFlags(data.globalStatusManager());

        double logA = M.sumOfLogsOfSquaredAmplitudes(data, parts);
        M.setParameterFlagsToUnchanged();

        LOG(INFO) << "logA = " << logA;

        std::cout << "Variable status after calculating:    \n";
        M.printFlags(data.globalStatusManager());

        /*if (gRandom->Uniform()>0.5) {
            double logATest = M.sumOfLogsOfSquaredAmplitudes(dataTest, partsTest);
            LOG(INFO) << "logATest = " << logATest;
            assert (logA == logATest);

            std::cout << "Variable status after calculating test:    \n";
            M.printFlags(data.globalStatusManager());

        }*/



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
