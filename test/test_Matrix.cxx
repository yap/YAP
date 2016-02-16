#include <catch.hpp>

#include <Constants.h>
#include <logging.h>
#include <LorentzTransformation.h>
#include <Matrix.h>
#include <Vector.h>

#include <cmath>

TEST_CASE( "Matrix" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);

    yap::ThreeMatrix<double> zero({0, 0, 0, 0, 0, 0, 0, 0, 0});
    yap::ThreeMatrix<double> unit({1, 0, 0, 0, 1, 0, 0, 0, 1});

    SECTION( "Initialization" ) {
        yap::ThreeMatrix<double> m = yap::zeroMatrix<double, 3>();
        REQUIRE(m == zero);

        auto u = yap::unitMatrix<double, 3>();
        REQUIRE(u == unit);
    }

    SECTION( "Transpose" ) {
        yap::ThreeMatrix<double> m({1, 2, 3, 4, 5, 6, 7, 8, 9});
        yap::ThreeMatrix<double> m_T({1, 4, 7, 2, 5, 8, 3, 6, 9});

        REQUIRE(m == m_T);
    }

    SECTION( "+ -" ) {
        yap::ThreeMatrix<double> m({1, 2, 3, 4, 5, 6, 7, 8, 9});
        yap::ThreeMatrix<double> minus_m({ -1, -2, -3, -4, -5, -6, -7, -8, -9});

        REQUIRE(m - m == zero);
        REQUIRE(-m == minus_m);
        REQUIRE(-1. * m == minus_m);
        REQUIRE(m + minus_m == zero);
        REQUIRE(m + m == 2. * m);
    }

    SECTION("boost") {
        yap::FourVector<double> a({1.2, 0., -1.1, 2.5});
        yap::FourVector<double> b({1.6, -0.2, 1.15, 1.5});
        yap::FourVector<double> c({1.9, 0.55, 2., -2.3});

        std::vector<yap::FourVector<double> > V({a, b, c});

        auto V_boost = -lorentzTransformation(V) * V;

        // check if they have been boosted into their rest frame
        REQUIRE( abs(std::accumulate(V_boost.begin(), V_boost.end(), yap::ThreeVector_0,
        [&](const yap::ThreeVector<double>& val, yap::FourVector<double> v) { return val + vect(v); })) == Approx(0.));

        // check if the invariant masses are the same
        for (unsigned i = 0; i < V.size(); ++i) {
            REQUIRE( norm(V[i]) == Approx(norm(V_boost[i])) );
        }
    }

}
