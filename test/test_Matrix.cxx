#include <catch.hpp>

#include <Exceptions.h>
#include <logging.h>
#include <LorentzTransformation.h>
#include <Matrix.h>
#include <Rotation.h>
#include <Vector.h>

#include <cmath>

TEST_CASE( "Matrix" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    const yap::ThreeMatrix<double> zero({0, 0, 0, 0, 0, 0, 0, 0, 0});
    const yap::ThreeMatrix<double> unit({1, 0, 0, 0, 1, 0, 0, 0, 1});

    SECTION( "Initialization" ) {
        const yap::ThreeMatrix<double> m = yap::zeroMatrix<double, 3>();
        REQUIRE(m == zero);

        auto u = yap::unitMatrix<double, 3>();
        REQUIRE(u == unit);

        const yap::ThreeMatrix<double> cov({1, 4, 6,  4, 2, 5,  6, 5, 3});
        const auto cov2 = yap::symmetricMatrix<double, 3>({1, 2, 3,  4, 5,  6});
        REQUIRE(cov == cov2);

        REQUIRE_THROWS_AS( (yap::symmetricMatrix<double, 2>({0, 1, 2, 3, 4})), yap::exceptions::Exception );
    }

    SECTION( "Transpose" ) {
        const yap::ThreeMatrix<double> m({1, 2, 3, 4, 5, 6, 7, 8, 9});
        const yap::ThreeMatrix<double> m_T({1, 4, 7, 2, 5, 8, 3, 6, 9});

        REQUIRE(transpose(m) == m_T);
    }

    SECTION( "+ -" ) {
        const yap::ThreeMatrix<double> m({1, 2, 3, 4, 5, 6, 7, 8, 9});
        const yap::ThreeMatrix<double> minus_m({ -1, -2, -3, -4, -5, -6, -7, -8, -9});

        REQUIRE(m - m == zero);
        REQUIRE(-m == minus_m);
        REQUIRE(-1. * m == minus_m);
        REQUIRE(m + minus_m == zero);
        REQUIRE(m + m == 2. * m);
    }

    SECTION( "Multiplication" ) {
        const yap::ThreeMatrix<double>  m({1, 2, 3, 4, 5, 6, 7, 8, 9});
        const yap::Matrix<double, 3, 4> n({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
        const yap::ThreeMatrix<double>  mm({30, 36, 42, 66, 81, 96, 102, 126, 150});
        const yap::Matrix<double, 3, 4> mn({38, 44, 50, 56, 83, 98, 113, 128, 128, 152, 176, 200});

        REQUIRE( unit * m == m );
        REQUIRE( -unit * m == -m );
        REQUIRE( m * unit == m );
        REQUIRE( (m * m) == mm );
        REQUIRE( (m * n) == mn );

        const auto v3 = yap::ThreeVector<double>({4, 3, 2});
        REQUIRE( m * v3 == yap::ThreeVector<double>({16, 43, 70}) );

        const yap::FourMatrix<double>  m4({1, 2, 3, 4,  5, 6, 7, 8,  9, 10, 11, 12,  13, 14, 15, 16});
        const auto v4 = yap::FourVector<double>({4, 3, 2, 1});
        REQUIRE( m * yap::FourVector<double>({4, 4, 3, 2}) == yap::FourVector<double>({4, 16, 43, 70}) );
        REQUIRE( m4 * v4 == yap::FourVector<double>({20, 60, 100, 140}) );

        // non-trivial template arguments
        yap::ThreeMatrix<yap::SquareMatrix<double, 2>> meh;
        meh[0][0] = yap::symmetricMatrix<double, 2>({1,2,3});
        meh[1][1] = yap::symmetricMatrix<double, 2>({4,5,6});
        meh[2][2] = yap::symmetricMatrix<double, 2>({7,8,9});

        REQUIRE( meh * unit == meh );
        REQUIRE( unit * meh == meh );
        REQUIRE( -unit * meh == -meh );
        REQUIRE( zero * meh == meh * zero );
        REQUIRE( meh - meh == meh * zero );
        REQUIRE( meh + meh == 2. * meh );
        REQUIRE( zero * meh != meh );
    }

    SECTION("boost") {
        /* //                              | this is the time component!!
         * // Example from ROOT            v
         * TLorentzVector a(0., -1.1, 2.5, 6.2);
         * a.Boost(-a.BoostVector());
         * a.Print()
         * (x,y,z,t)=(0.000000,-0.000000,0.000000,5.565968) (P,eta,phi,E)=(0.000000,1.443635,-1.570796,5.565968)
         */

        const yap::FourVector<double> a({6.2,  0.,  -1.1,  2.5});
        const yap::FourVector<double> b({8.6, -0.2,  1.15, 1.5});
        const yap::FourVector<double> c({6.9,  0.55, 2.,  -2.3});

        // see if we can boost a 4vector into its rest frame
        REQUIRE( abs(vect(lorentzTransformation(boost(-a)) * a)) == Approx(0.) );
        REQUIRE( norm(lorentzTransformation(boost(-a)) * a) == Approx(norm(a)) );

        REQUIRE( abs(vect(lorentzTransformation(-boost(a)) * a)) == Approx(0.) );
        REQUIRE( norm(lorentzTransformation(-boost(a)) * a) == Approx(norm(a)) );

        REQUIRE( abs(vect(lorentzTransformation(-a) * a)) == Approx(0.) );
        REQUIRE( norm(lorentzTransformation(-a) * a) == Approx(norm(a)) );

        // now boost all 4vectors
        std::vector<yap::FourVector<double> > V({a, b, c});
        const auto V_boosted = lorentzTransformation(-V) * V;

        // check if they have been boosted into their rest frame
        REQUIRE( abs(vect(std::accumulate(V_boosted.begin(), V_boosted.end(), yap::FourVector<double>()))) == Approx(0.));

        // check if the invariant masses are the same
        for (unsigned i = 0; i < V.size(); ++i) {
            REQUIRE( norm(V[i]) == Approx(norm(V_boosted[i])) );
        }
    }

    SECTION("ThreeVector rotations") {

        for (double alpha = 0; alpha < 3.14; alpha += 0.2) {

            for (unsigned i = 0; i < 3; ++i) {
                yap::ThreeVector<double> axis;
                axis[i] = 1;

                const yap::ThreeVector<double> a({0., -1.1, 2.5});
                const yap::ThreeVector<double> b({ -0.2, 1.15, 1.5});
                const yap::ThreeVector<double> c({0.55, 2., -2.3});

                const yap::ThreeVector<double> x({1., 0., 0.});
                const yap::ThreeVector<double> y({0., 1., 0.});
                const yap::ThreeVector<double> z({0., 0., 1.});

                const auto trans = yap::rotation<double>(axis, alpha);

                const auto a_trans = trans * a;
                const auto b_trans = trans * b;
                const auto c_trans = trans * c;

                REQUIRE( norm(a_trans) == Approx(norm(a)) );
                REQUIRE( norm(b_trans) == Approx(norm(b)) );
                REQUIRE( norm(c_trans) == Approx(norm(c)) );

                REQUIRE( angle((a_trans), (b_trans)) == Approx(angle((a), (b))) );
                REQUIRE( angle((a_trans), (c_trans)) == Approx(angle((a), (c))) );
                REQUIRE( angle((b_trans), (c_trans)) == Approx(angle((b), (c))) );

                switch (i) {
                    case 0:
                        REQUIRE( angle((trans * x), (x)) == Approx(0.) );
                        REQUIRE( angle((trans * y), (y)) == Approx(alpha) );
                        REQUIRE( angle((trans * z), (z)) == Approx(alpha) );
                        break;

                    case 1:
                        REQUIRE( angle((trans * x), (x)) == Approx(alpha) );
                        REQUIRE( angle((trans * y), (y)) == Approx(0.) );
                        REQUIRE( angle((trans * z), (z)) == Approx(alpha) );
                        break;

                    case 2:
                        REQUIRE( angle((trans * x), (x)) == Approx(alpha) );
                        REQUIRE( angle((trans * y), (y)) == Approx(alpha) );
                        REQUIRE( angle((trans * z), (z)) == Approx(0.) );
                        break;
                }


                const auto zeroTrans = yap::rotation<double>(axis, 0);

                REQUIRE ( a == zeroTrans * a );
                REQUIRE ( b == zeroTrans * b );
                REQUIRE ( c == zeroTrans * c );

            }
        }

    }

    SECTION("FourVector rotations") {

        for (double alpha = 0; alpha < 3.14; alpha += 0.2) {

            for (unsigned i = 0; i <3; ++i) {

                yap::ThreeVector<double> axis;
                axis[i] = 1;

                const yap::FourVector<double> a({6.2, 0., -1.1, 2.5});
                const yap::FourVector<double> b({8.6, -0.2, 1.15, 1.5});
                const yap::FourVector<double> c({6.9, 0.55, 2., -2.3});

                const yap::FourVector<double> x({6.2, 1., 0., 0.});
                const yap::FourVector<double> y({6.2, 0., 1., 0.});
                const yap::FourVector<double> z({6.9, 0., 0., 1.});

                const auto trans = lorentzTransformation( yap::rotation<double>(axis, alpha) );

                const auto a_trans = trans * a;
                const auto b_trans = trans * b;
                const auto c_trans = trans * c;

                REQUIRE( norm(a_trans) == Approx(norm(a)) );
                REQUIRE( norm(b_trans) == Approx(norm(b)) );
                REQUIRE( norm(c_trans) == Approx(norm(c)) );

                REQUIRE( angle(vect(a_trans), vect(b_trans)) == Approx(angle(vect(a), vect(b))) );
                REQUIRE( angle(vect(a_trans), vect(c_trans)) == Approx(angle(vect(a), vect(c))) );
                REQUIRE( angle(vect(b_trans), vect(c_trans)) == Approx(angle(vect(b), vect(c))) );

                switch (i) {
                    case 0:
                        REQUIRE( angle(vect(trans * x), vect(x)) == Approx(0.) );
                        REQUIRE( angle(vect(trans * y), vect(y)) == Approx(alpha) );
                        REQUIRE( angle(vect(trans * z), vect(z)) == Approx(alpha) );
                        break;

                    case 1:
                        REQUIRE( angle(vect(trans * x), vect(x)) == Approx(alpha) );
                        REQUIRE( angle(vect(trans * y), vect(y)) == Approx(0.) );
                        REQUIRE( angle(vect(trans * z), vect(z)) == Approx(alpha) );
                        break;

                    case 2:
                        REQUIRE( angle(vect(trans * x), vect(x)) == Approx(alpha) );
                        REQUIRE( angle(vect(trans * y), vect(y)) == Approx(alpha) );
                        REQUIRE( angle(vect(trans * z), vect(z)) == Approx(0.) );
                        break;
                }


                const auto zeroTrans = lorentzTransformation( yap::rotation<double>(axis, 0) );

                REQUIRE ( a == zeroTrans * a );
                REQUIRE ( b == zeroTrans * b );
                REQUIRE ( c == zeroTrans * c );

            }
        }

    }
}
