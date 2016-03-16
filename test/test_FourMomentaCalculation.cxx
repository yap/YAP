#include <catch.hpp>
#include <catch_capprox.hpp>

#include <Constants.h>
#include <FinalStateParticle.h>
#include <DecayingParticle.h>
#include <FourVector.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <MassAxes.h>
#include <make_unique.h>
#include <Model.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <ZemachFormalism.h>

#include <cmath>

bool valid_5d(const double m2_12, const double m2_14, const double m2_23,
        const double m2_34, const double m2_13,
        const double m_Parent,
        const double m_a, const double m_b,
        const double m_c, const double m_d)
{
    // 5D hypercube lower boundaries
    if ( m2_12 < pow(m_a + m_b, 2) ||
            m2_14 < pow(m_a + m_d, 2) ||
            m2_23 < pow(m_b + m_c, 2) ||
            m2_34 < pow(m_c + m_d, 2) ||
            m2_13 < pow(m_a + m_c, 2) ) {
        return false;
    }

    // 5D hypercube upper boundaries
    if ( m2_12 > pow(m_Parent - m_c - m_d, 2) ||
            m2_14 > pow(m_Parent - m_b - m_c, 2) ||
            m2_23 > pow(m_Parent - m_a - m_d, 2) ||
            m2_34 > pow(m_Parent - m_a - m_b, 2) ||
            m2_13 > pow(m_Parent - m_b - m_d, 2) ) {
        return false;
    }

    const double m2_24 = (m_Parent*m_Parent + 2.*(m_a*m_a + m_b*m_b + m_c*m_c + m_d*m_d)) - m2_12 - m2_14 - m2_23 - m2_34 - m2_13;

    if ( sqrt(m2_12) + sqrt(m2_34) > m_Parent ||
            sqrt(m2_14) + sqrt(m2_23) > m_Parent ||
            sqrt(m2_13) + sqrt(m2_24) > m_Parent ) {
        return false;
    }

    const double M = m2_12;
    const double M2 = M*M;

    const double N = m2_34;
    const double N2 = N*N;

    const double P = m2_12 + m2_14 + m2_24 - m_a*m_a - m_b*m_b - m_d*m_d; // m2_124
    const double P2 = P*P;

    const double Q = m2_13 + m2_14 + m2_34 - m_a*m_a - m_c*m_c - m_d*m_d; // m2_134;
    const double Q2 = Q*Q;

    const double R = m2_14;
    const double R2 = R*R;

    const double m = m_c*m_c;
    const double m2 = m*m;

    const double n = m_b*m_b;
    const double n2 = n*n;

    const double p = m_a*m_a;
    const double p2 = p*p;

    const double q = m_d*m_d;
    const double q2 = q*q;

    const double r = m_Parent*m_Parent; // E^2
    const double r2 = r*r;

    double B = (M2*Q2 + N2*P2 + M2*R2 + N2*R2 + P2*Q2) - 2.*(M2*Q*R + N2*P*R + M*N*R2 + M*P*Q2 + N*P2*Q)
                + 2.*(M*N*P*Q + M*N*P*R + M*N*Q*R + M*P*Q*R + N*P*Q*R) - 2.*(M2*Q*m + N2*P*n + M2*R*m + N2*R*n + M*Q2*q
                        + N*P2*p + M*R2*r + N*R2*r + P2*Q*p + P*Q2*q) - 2.*(M*N*P*m + M*N*Q*n + M*P*R*p + N*Q*R*q + P*Q*R*r)
                        + 2.*(M*N*R*(m + n - 2.*r) + M*P*Q*(m + p - 2.*q) + N*Q*P*(n + q - 2.*p) + Q*R*M*(q + r - 2.*m) + P*R*N*(p + r - 2.*n))
                        + (M2*m2 + N2*n2 + P2*p2 + Q2*q2 + R2*r2) + 2.*(M*N*m*n + M*P*m*p + N*Q*n*q + P*R*p*r + Q*R*q*r)
                        + 2.*(M*Q*(m*q + m*n + q*n + m*p + q*r - p*r) + N*P*(n*p + n*m + p*m + p*r + n*q - q*r) + M*R*(m*r + m*p + r*p + m*n + r*q - n*q)
                                + N*R*(n*r + n*q + r*q + n*m + r*p - m*p) + P*Q*(p*q + p*r + q*r + p*m + q*n - m*n) ) - 2.*(M*m*(m*p + m*n + q*r - p*r - n*q + 2.*n*p)
                                        + N*n*(n*m + n*q + p*r - p*m - q*r + 2.*m*q) + P*p*(p*m + p*r + n*q - m*n - q*r + 2.*m*r) + Q*q*(q*n + q*r + m*p - m*n - p*r + 2.*n*r)
                                        + R*r*(r*p + r*q + m*n - m*p - n*q + 2.*p*q) ) + (m2*n2 + m2*p2 + n2*q2 + p2*r2 + q2*r2) - 2.*(m2*n*p + m*n2*q + m*p2*r + n*q2*r + p*q*r2)
                                        + 2.*(m*n*p*q + m*n*p*r + m*n*q*r + m*p*q*r + n*p*q*r);

    if (B < 0.)
        return true;

    return false;
}


TEST_CASE( "FourMomentaCalculation" )
{

    SECTION("3 final state particles") {


        // disable logs in text
        yap::disableLogs(el::Level::Global);

        // load particle factory
        yap::ParticleFactory factory("../data/evt.pdl");

        // create final state particles
        auto piPlus = factory.fsp(211);
        auto KMinus = factory.fsp(-321);

        // create model and set initial and final state
        yap::Model M(std::make_unique<yap::ZemachFormalism>());
        M.setFinalState({piPlus, KMinus, piPlus});

        // create book-keeping resonance
        /// \todo Allow direct phase-space three-body decay
        auto X = factory.decayingParticle(factory.pdgCode("f_0"), 3);
        X->addChannel({piPlus, KMinus});

        // create initial state particle
        auto D = factory.decayingParticle(factory.pdgCode("D+"), 3);
        D->addChannel({X, piPlus});

        // REQUIRE( M.consistent() );

        // choose Dalitz coordinates m^2_12 and m^2_23
        const yap::MassAxes massAxes = M.massAxes({{0, 1}, {1, 2}});

        auto m_0_range = M.massRange(massAxes[0]);
        auto m_1_range = M.massRange(massAxes[1]);

        const unsigned N = 200;

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
                        if (std::find(massAxes[1]->indices().begin(), massAxes[1]->indices().end(), i) == massAxes[1]->indices().end())
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


    SECTION("4 final state particles") {

        // disable logs in text
        yap::disableLogs(el::Level::Global);

        // load particle factory
        yap::ParticleFactory factory("../data/evt.pdl");

        // create final state particles
        auto piPlus = factory.fsp(211);
        auto piMinus = factory.fsp(-211);
        auto KPlus = factory.fsp(321);
        auto KMinus = factory.fsp(-321);

        // create model and set initial and final state
        yap::Model M(std::make_unique<yap::HelicityFormalism>());
        M.setFinalState({piPlus, piMinus, KPlus, KMinus});

        // create book-keeping resonance
        /// \todo Allow direct phase-space three-body decay
        auto X = factory.decayingParticle(factory.pdgCode("f_0"), 3);
        X->addChannel({piPlus, piMinus});

        auto X2 = factory.decayingParticle(factory.pdgCode("f_0"), 3);
        X2->addChannel({KPlus, KMinus});

        // create initial state particle
        auto D = factory.decayingParticle(factory.pdgCode("D0"), 3);
        D->addChannel({X, X2});

        // REQUIRE( M.consistent() );

        // choose Dalitz coordinates m^2_12, m^2_14, m^2_23, m^2_34, m^2_13
        const yap::MassAxes massAxes = M.massAxes({{0, 1}, {0, 3}, {1, 2}, {2, 3}, {0, 2}});

        auto m_0_range = M.massRange(massAxes[0]);
        auto m_1_range = M.massRange(massAxes[1]);
        auto m_2_range = M.massRange(massAxes[2]);
        auto m_3_range = M.massRange(massAxes[3]);
        auto m_4_range = M.massRange(massAxes[4]);

        const unsigned N = 15;

        for (double m_0 = m_0_range[0]; m_0 <= m_0_range[1]; m_0 += (m_0_range[1] - m_0_range[0]) / N) {
            for (double m_1 = m_1_range[0]; m_1 <= m_1_range[1]; m_1 += (m_1_range[1] - m_1_range[0]) / N)
                for (double m_2 = m_2_range[0]; m_2 <= m_2_range[1]; m_2 += (m_2_range[1] - m_2_range[0]) / N)
                    for (double m_3 = m_3_range[0]; m_3 <= m_3_range[1]; m_3 += (m_3_range[1] - m_3_range[0]) / N)
                        for (double m_4 = m_4_range[0]; m_4 <= m_4_range[1]; m_4 += (m_4_range[1] - m_4_range[0]) / N) {

                            std::vector<double> m2 = {m_0 * m_0, m_1 * m_1, m_2 * m_2, m_3 * m_3, m_4 * m_4};

                            // calculate final state momenta
                            auto P = M.calculateFourMomenta(massAxes, m2);

                            //-------------------------
                            // check phase space
                            bool inPhaseSpace = valid_5d(m2[0], m2[1], m2[2], m2[3], m2[4],
                                    M.initialStateParticle()->mass()->value(),
                                    M.finalStateParticles()[0]->mass()->value(),
                                    M.finalStateParticles()[1]->mass()->value(),
                                    M.finalStateParticles()[2]->mass()->value(),
                                    M.finalStateParticles()[3]->mass()->value());

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
                                REQUIRE( abs(p_isp) == Approx(M.initialStateParticle()->mass()->value()) );

                                // check Dalitz axes
                                for (size_t i = 0; i < massAxes.size(); ++i) {
                                    auto p = std::accumulate(massAxes[i]->indices().begin(), massAxes[i]->indices().end(), yap::FourVector_0, [&](const yap::FourVector<double>& p, unsigned j) {return p + P[j];});
                                    REQUIRE( norm(p) == Approx(m2[i]) );
                                }

                            }
                        }
        }
    }
}

