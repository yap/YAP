#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <Exceptions.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <MathUtilities.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <Resonance.h>
#include <WignerD.h>
#include <ZemachFormalism.h>

#include <TGenPhaseSpace.h>
#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <cmath>

/*
 * Test the calculation of helicity angles
 */

// copied from rootPWA
TLorentzRotation hfTransform(const TLorentzVector& daughterLv)
{
    TLorentzVector daughter = daughterLv;
    const TVector3 zAxisParent(0, 0, 1);  // take z-axis as defined in parent frame
    const TVector3 yHfAxis = zAxisParent.Cross(daughter.Vect());  // y-axis of helicity frame
    // rotate so that yHfAxis becomes parallel to y-axis and zHfAxis ends up in (x, z)-plane
    TRotation rot1;
    rot1.RotateZ(0.5*yap::pi<double>() - yHfAxis.Phi());
    rot1.RotateX(yHfAxis.Theta() - yap::pi<double>());
    daughter *= rot1;
    // rotate about yHfAxis so that daughter momentum is along z-axis
    TRotation rot2;
    rot2.RotateY(-yap::signum(daughter.X()) * daughter.Theta());
    daughter *= rot2;
    // boost to daughter RF
    rot1.Transform(rot2);
    TLorentzRotation hfTransform(rot1);
    hfTransform.Boost(-daughter.BoostVector());
    return hfTransform;
}

void transformDaughters(const yap::Model& M,
        std::map<const std::shared_ptr<yap::ParticleCombination>, std::vector<double> >& results,
        const std::shared_ptr<yap::ParticleCombination>& pc,
        std::vector<TLorentzVector> finalStatesHf)
{
    // loop over daughters
    for (const auto& daugh : pc->daughters()) {

        if (daugh->daughters().empty())
            continue;

        // construct 4-vector of daughter
        TLorentzVector daughter;
        for (unsigned i : daugh->indices())
            daughter += finalStatesHf.at(i);

        results[daugh] = {daughter.Phi(), daughter.Theta()};

        // next helicity frame
        const TLorentzRotation transDaugh = hfTransform(daughter);
        for (unsigned i : daugh->indices())
            finalStatesHf.at(i).Transform(transDaugh);

        transformDaughters(M, results, daugh, finalStatesHf);
    }
}


TEST_CASE( "HelicityAngles" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);
    //yap::plainLogs(el::Level::Debug);

    // init random generator
    gRandom->SetSeed(1234);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    yap::Model M(std::make_unique<yap::ZemachFormalism>());

    yap::ParticleFactory factory((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    // initial state particle
    auto D = factory.decayingParticle(factory.pdgCode("D+"), radialSize);

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // set final state
    M.setFinalState({piPlus, piMinus, piPlus});

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::BreitWigner>());
    rho->addChannel({piPlus, piMinus});

    // Add channels to D
    D->addChannel({rho,      piPlus});

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = M.getMassAxes({{0, 1}, {1, 2}});

    REQUIRE( M.consistent() );

    const int nParticles(3);

    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    Double_t masses[nParticles] = { piPlus->mass()->value(), piMinus->mass()->value(),
                                    piPlus->mass()->value()
                                  };

    for (unsigned int iEvt = 0; iEvt < 1; ++iEvt) {
        TGenPhaseSpace event;
        event.SetDecay(P, nParticles, masses);
        event.Generate();

        std::vector<yap::FourVector<double> > momenta;
        std::vector<TLorentzVector> rootMomenta;

        for (unsigned i = 0; i < nParticles; ++i) {
            TLorentzVector p = *event.GetDecay(i);
            rootMomenta.push_back(p);
            momenta.push_back(yap::FourVector<double>({p.T(), p.X(), p.Y(), p.Z()}));
        }

        M.addDataPoint(momenta);
        auto dp = M.dataSet().back();

        // calculate helicity angles with alternative method copied from rootPWA
        std::map<const std::shared_ptr<yap::ParticleCombination>, std::vector<double> > results;

        for (auto pc : D->particleCombinations()) {
            REQUIRE( M.fourMomenta()->m(dp, pc) == Approx(D->mass()->value()) );

            transformDaughters(M, results, pc, rootMomenta);
        }

        // compare results
        for (auto& kv : results) {
            REQUIRE( M.helicityAngles()->phi(dp, kv.first) == Approx(results.at(kv.first)[0]) );
            REQUIRE( M.helicityAngles()->theta(dp, kv.first) == Approx(results.at(kv.first)[1]) );
        }


    }
}
