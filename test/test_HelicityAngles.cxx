#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <LorentzTransformation.h>
#include <MassAxes.h>
#include <Model.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <Resonance.h>

#include <TGenPhaseSpace.h>
#include <TLorentzRotation.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <assert.h>
#include <cmath>

/*
 * Test the calculation of helicity angles
 */

// copied from rootPWA
TLorentzRotation hfTransform(const TLorentzVector& daughterLv)
{
    std::cout << "original daughter: "; daughterLv.Print();
    TLorentzVector daughter = daughterLv;
    const TVector3 zAxisParent(0, 0, 1);  // take z-axis as defined in parent frame
    const TVector3 yHfAxis = zAxisParent.Cross(daughter.Vect());  // y-axis of helicity frame
    // rotate so that yHfAxis becomes parallel to y-axis and zHfAxis ends up in (x, z)-plane
    TRotation rot1;
    rot1.RotateZ(0.5*yap::pi<double>() - yHfAxis.Phi());
    //rot1.RotateX(yHfAxis.Theta() - yap::pi<double>() / 2.);
    daughter *= rot1;
    // rotate about yHfAxis so that daughter momentum is along z-axis
    TRotation rot2;
    rot2.RotateY(-yap::signum(daughter.X()) * daughter.Theta());
    daughter *= rot2;
    std::cout << "rotated daughter: "; daughter.Print();
    // boost to daughter RF
    rot1.Transform(rot2);
    TLorentzRotation hfTransform(rot1);
    hfTransform.Boost(-daughter.BoostVector());
    return hfTransform;
}

void transformDaughters(const yap::Model& M,
        std::map<const std::shared_ptr<yap::ParticleCombination>, std::array<double, 2> >& results,
        const std::shared_ptr<yap::ParticleCombination>& pc,
        std::vector<TLorentzVector> finalStatesHf)
{
    // loop over daughters
    for (const auto& daugh : pc->daughters()) {

        //if (daugh->daughters().empty())
          //  continue;

        // construct 4-vector of daughter
        TLorentzVector daughter;
        for (unsigned i : daugh->indices())
            daughter += finalStatesHf.at(i);

        results[pc] = {daughter.Phi(), daughter.Theta()};

        // next helicity frame
        const TLorentzRotation transDaugh = hfTransform(daughter);
        for (unsigned i : daugh->indices())
            finalStatesHf.at(i).Transform(transDaugh);

        transformDaughters(M, results, daugh, finalStatesHf);
    }
}


// YAP version
yap::FourMatrix<double> transformation_to_helicityFrame(const yap::FourVector<double>& daughter)
{
    // rotate to put x-y component of daughter in y direction
    // rotate to put daughter in z direction
    auto R = rotation(yap::ThreeAxis_Y, -yap::theta(vect(daughter), yap::ThreeAxes))
            * rotation(yap::ThreeAxis_Z, -yap::phi(vect(daughter), yap::ThreeAxes));

    // return boost * rotations
    return lorentzTransformation(-(R * daughter)) * lorentzTransformation(R);
}

void calculate_helicity_angles(const yap::Model& M,
                               std::map<const std::shared_ptr<yap::ParticleCombination>, std::array<double, 2> >& phi_theta,
                               const std::shared_ptr<yap::ParticleCombination>& pc,
                               std::vector<yap::FourVector<double> > momenta)
{
    // loop over daughters
    for (const auto& d : pc->daughters()) {

        auto p = yap::FourVector_0;
        // construct 4-vector of daughter
        for (const auto& i : d->indices())
            p += momenta[i];

        // if not yet calculated
        if (phi_theta.find(pc) == phi_theta.end()) {
            std::cout << "calculate angles for daughter " << to_string(*d) << ": " << to_string(p) <<"\n";
            phi_theta[pc] = angles(vect(p), yap::ThreeAxes);

            std::cout << "(phi, theta) = (" << phi_theta[pc][0] << ", " << phi_theta[pc][1] << ")\n";
        }

        /*else {
            // check that results would be the same within numerical uncertainty
            for (unsigned i = 0; i < 2; ++i) {
                std::cout << "i " << i << " new " << hAngles[i] << "; old " << phi_theta[pc][i] << std::endl;
                std::cout << fabs(fabs(hAngles[i] - phi_theta[pc][i]) - yap::pi<double>()) << std::endl;
                std::cout << fabs(fabs(hAngles[i] + phi_theta[pc][i]) - yap::pi<double>()) << std::endl;
                assert((i == 0 && fabs(fabs(hAngles[i] - phi_theta[pc][i]) - yap::pi<double>()) < 1e-10) ||
                       (i == 1 && fabs(fabs(hAngles[i] + phi_theta[pc][i]) - yap::pi<double>()) < 1e-10) );
            }
        }*/

        if (d->daughters().empty())
            continue;

        // next helicity frame
        const auto L = transformation_to_helicityFrame(p);

        std::cout << "transforming for    " << to_string(*d) << "\n";
        std::cout << "with transformation " << to_string(L) << "\n";
        for (unsigned i : d->indices()) {
            std::cout << i << " before: " << to_string(momenta[i]) << "\n";
            momenta[i] = L * momenta[i];
            std::cout << i << " after:  " << to_string(momenta[i]) << "\n";
        }

        calculate_helicity_angles(M, phi_theta, d, momenta);
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

    yap::Model M(std::make_unique<yap::HelicityFormalism>());

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
    const yap::MassAxes massAxes = M.massAxes({{0, 1}, {1, 2}});

    REQUIRE( M.consistent() );

    // create DataSet
    auto data = M.dataSet();

    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    std::vector<double> masses = { piPlus->mass()->value(), piMinus->mass()->value(), piPlus->mass()->value() };

    for (unsigned int iEvt = 0; iEvt < 100; ++iEvt) {
        TGenPhaseSpace event;
        event.SetDecay(P, masses.size(), &masses[0]);
        event.Generate();

        std::vector<yap::FourVector<double> > momenta;
        std::vector<TLorentzVector> rootMomenta;

        for (unsigned i = 0; i < masses.size(); ++i) {
            TLorentzVector p = *event.GetDecay(i);

            // just testing. Theta of downstream helicity angles must stay the same
            // Phi can change, since it only affects the phase of the amplitude, and it will change in the same way for all amplitudes

            p.Boost(0.1, 0., 0.); // changes YAP and ROOT in the same way, downstream angles are the SAME as unboosted. M is wrong.
            //p.Boost(-0.1, 0., 0.);

            //p.RotateX(0.25); // changes YAP and ROOT in the same way, downstream angles are NOT the same as unrotated. M is wrong
            //p.RotateY(0.25); // changes YAP and ROOT in the same way, downstream angles are NOT the same as unrotated. M is wrong

            //p.RotateZ(0.25); // changes M, YAP and ROOT in the same way, downstream angles are the SAME as unrotated.

            rootMomenta.push_back(p);
            momenta.push_back(yap::FourVector<double>({p.T(), p.X(), p.Y(), p.Z()}));
        }

        TLorentzVector p;
        for (auto& m : rootMomenta) {
            p += m;
        }
        std::cout << "ISP momentum = "; p.Print();


        data.add(momenta);
        const auto dp = data.points().back();

        std::map<const std::shared_ptr<yap::ParticleCombination>, std::array<double, 2> > phi_theta;
        std::map<const std::shared_ptr<yap::ParticleCombination>, std::array<double, 2> > phi_thetaRoot;

        for (auto pc : D->particleCombinations()) {
            REQUIRE( M.fourMomenta()->m(dp, pc) == Approx(D->mass()->value()) );

            /// \todo they do not give the same results
            calculate_helicity_angles(M, phi_theta, pc, momenta); // YAP
            transformDaughters(M, phi_thetaRoot, pc, rootMomenta); // rootPWA
        }

        // compare results
        for (auto& kv : phi_theta) {
            std::cout << yap::to_string(*kv.first) << "\n";
            std::cout << "M.helicityAngles():  (" <<  M.helicityAngles()->phi(dp, kv.first) << ", " << M.helicityAngles()->theta(dp, kv.first) << ")\n";
            std::cout << "YAP  helicityAngles: (" <<  kv.second[0] << ", " << kv.second[1] << ")\n";
            std::cout << "root helicityAngles: (" <<  phi_thetaRoot[kv.first][0] << ", " << phi_thetaRoot[kv.first][1] <<
                    ") equiv (" << ((phi_thetaRoot[kv.first][0] > yap::pi<double>()/2.) ? (phi_thetaRoot[kv.first][0] - yap::pi<double>()) : (phi_thetaRoot[kv.first][0] + yap::pi<double>())) <<
                            ", " <<  yap::pi<double>() - phi_thetaRoot[kv.first][1] << ")\n";

            //REQUIRE( M.helicityAngles()->phi(dp, kv.first)   == Approx(kv.second[0]) );
            //REQUIRE( M.helicityAngles()->theta(dp, kv.first) == Approx(kv.second[1]) );
        }


        std::cout << "\n";


    }
}
