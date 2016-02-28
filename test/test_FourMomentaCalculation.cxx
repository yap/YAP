#include <catch.hpp>
#include <catch_capprox.hpp>

#include <Constants.h>
#include <FinalStateParticle.h>
#include <DecayingParticle.h>
#include <FourVector.h>
#include <logging.h>
#include <MassAxes.h>
#include <make_unique.h>
#include <Model.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <ZemachFormalism.h>

#include <cmath>

TEST_CASE( "FourMomentaCalculation" )
{

    // disable logs in text
    // yap::disableLogs(el::Level::Global);

    // load particle factory
    yap::ParticleFactory factory("../data/evt.pdl");

    // create final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // create model and set initial and final state
    yap::Model M(std::make_unique<yap::ZemachFormalism>());
    M.setFinalState({piPlus, piMinus, piPlus});

    // create book-keeping resonance
    /// \todo Allow direct phase-space three-body decay
    auto X = factory.decayingParticle(factory.pdgCode("f_0"), 3);
    X->addChannel({piPlus, piMinus});

    // create initial state particle
    auto D = factory.decayingParticle(factory.pdgCode("D+"), 3);
    D->addChannel({X, piPlus});

    // REQUIRE( M.consistent() );

    // choose Dalitz coordinates m^2_12 and m^2_23
    const yap::MassAxes massAxes = M.massAxes({{0, 1}, {1, 2}});

    auto m_0_range = M.massRange(massAxes[0]);
    auto m_1_range = M.massRange(massAxes[1]);

    const unsigned N = 20;
    for (double m_0 = m_0_range[0]; m_0 <= m_0_range[1]; m_0 += (m_0_range[1] - m_0_range[0]) / N) {
        for (double m_1 = m_1_range[0]; m_1 <= m_1_range[1]; m_1 += (m_1_range[1] - m_1_range[0]) / N) {

            std::vector<double> m2 = {m_0 * m_0, m_1 * m_1};

            // calculate final state momenta
            auto P = M.calculateFourMomenta(massAxes, m2);

            //-------------------------
            // check phase space
            double m_isp = M.initialStateParticle()->mass()->value();
            // find a, b, and c such that mass axes are (ab) (bc)
            double m_a = 0;
            double m_b = 0;
            double m_c = 0;
            for (unsigned i = 0; i < M.finalStateParticles().size(); ++i) {
                if (std::find(massAxes[0]->indices().begin(), massAxes[0]->indices().end(), i) != massAxes[0]->indices().end()) {
                    if (std::find(massAxes[1]->indices().begin(), massAxes[1]->indices().end(), i) == massAxes[0]->indices().end())
                        m_a = M.finalStateParticles()[i]->mass()->value();
                    else
                        m_b = M.finalStateParticles()[i]->mass()->value();
                } else
                    m_c = M.finalStateParticles()[i]->mass()->value();
            }
            // check constraints on m2[0]
            bool inPhaseSpace = m2[0] >= pow(m_a + m_b, 2) and m2[0] <= pow(m_isp - m_c, 2);
            // check constraints on m2[1]
            if (inPhaseSpace) {
                double Eb = (m2[0] - m_a * m_a + m_b * m_b) / 2. / sqrt(m2[0]);
                double Ec = (m_isp * m_isp - m2[0] - m_c * m_c) / 2. / sqrt(m2[0]);
                double Pb = sqrt(Eb * Eb - m_b * m_b);
                double Pc = sqrt(Ec * Ec - m_c * m_c);
                inPhaseSpace = fabs(m2[1] - m_b * m_b - m_c * m_c - 2. * Eb * Ec) <= 2. * Pb * Pc;
            }

            // require P is empty if outside phase space
            REQUIRE( P.empty() == !inPhaseSpace );

            if (!P.empty()) {

                //-------------------------
                // check fsp masses
                for (size_t i = 0; i < P.size(); ++i)
                    REQUIRE( abs(P[i]) == Approx(M.finalStateParticles()[i]->mass()->value()) );

                //-------------------------
                // check isp mass
                // isp
                auto p_isp = std::accumulate(P.begin(), P.end(), yap::FourVector_0);
                REQUIRE( abs(p_isp) == Approx(m_isp) );

                // check Dalitz axes
                for (size_t i = 0; i < massAxes.size(); ++i) {
                    auto p = std::accumulate(massAxes[i]->indices().begin(), massAxes[i]->indices().end(), yap::FourVector_0, [&](const yap::FourVector<double>& p, unsigned j) {return p + P[j];});
                    REQUIRE( norm(p) == Approx(m2[i]) );
                }

            }
        }
    }
}
