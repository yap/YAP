#include <catch.hpp>
#include <catch_capprox.hpp>

#include <AmplitudeBasis.h>
#include <ComplexBasis.h>
#include <Constants.h>
#include <logging.h>
#include <Matrix.h>

template <size_t N>
inline bool compare_covariances(const yap::SquareMatrix<double, N>& c1, const yap::SquareMatrix<double, N>& c2)
{
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            if (not (c1[i][j] == Approx(c2[i][j])))
                return false;
    return true;
}

inline bool compare_covariances(const yap::SquareMatrix<yap::SquareMatrix<double, 2>, 3>& c1, const yap::SquareMatrix<yap::SquareMatrix<double, 2>, 3>& c2)
{
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            if (not compare_covariances(c1[i][j], c2[i][j]))
                return false;
    return true;
}

inline bool compare_complex_basis(const yap::complex_basis::basis<double>& c1, const yap::complex_basis::basis<double>& c2)
{
    if (not (std::complex<double>(c1) == Catch::Detail::CApprox(std::complex<double>(c2))))
        return false;

    if (not compare_covariances(c1.covariance(), c2.covariance()) )
        return false;

    return true;
}

TEST_CASE( "ComplexBasisTransformations" )
{
    SECTION( "const" ) {
        const yap::complex_basis::polar<double>     p1(2.1, yap::rad(73.), {0.1, yap::rad(10.), 0.01});

        const yap::complex_basis::cartesian<double> c1(p1);
        const yap::complex_basis::polar<double>     p2(c1);

        // check forth-back conversion
        REQUIRE( std::complex<double>(p1) == Catch::Detail::CApprox(std::complex<double>(p2)) );
        REQUIRE( p1.value() == p2.value() );
        REQUIRE( compare_covariances(p1.covariance(), p2.covariance()) );

        // check if value & cov are different in the two bases
        REQUIRE( not (p1.value() == c1.value()) );
        REQUIRE( not compare_covariances(p1.covariance(), c1.covariance()) );
        REQUIRE( p1.covariance()[0][0] != c1.covariance()[0][0]);

        // check casting
        REQUIRE( std::complex<double>(p1) == Catch::Detail::CApprox(std::complex<double>(c1)) );
    }
    SECTION( "non-const" ) {
        yap::complex_basis::polar<double>     p1(1.1, yap::rad(73.), {0.1, yap::rad(10.), 0.01});

        yap::complex_basis::cartesian<double> c1(p1);
        yap::complex_basis::polar<double>     p2(c1);

        // check forth-back conversion
        REQUIRE( std::complex<double>(p1) == Catch::Detail::CApprox(std::complex<double>(p2)) );
        REQUIRE( p1.value() == p2.value() );
        REQUIRE( compare_covariances(p1.covariance(), p2.covariance()) );

        // check if value & cov are different in the two bases
        REQUIRE( not (p1.value() == c1.value()) );
        REQUIRE( not compare_covariances(p1.covariance(), c1.covariance()) );
        REQUIRE( p1.covariance()[0][0] != c1.covariance()[0][0]);

        // check casting
        REQUIRE( std::complex<double>(p1) == Catch::Detail::CApprox(std::complex<double>(c1)) );
    }
}

TEST_CASE( "SpinBasisTransformations" )
{
    SECTION( "const" ) {
        // amplitudes in polar coords
        const yap::complex_basis::polar<double> pol_1(2.624, yap::rad( 47.), {1.210, 1.951, 0.101});
        const yap::complex_basis::polar<double> pol_2(3.157, yap::rad(120.), {2.176, 1.686, 0.153});
        const yap::complex_basis::polar<double> pol_3(4.384, yap::rad(163.), {1.782, 2.064, 0.279});

        const yap::complex_basis::cartesian<double> cart_1(pol_1);
        const yap::complex_basis::cartesian<double> cart_2(pol_2);
        const yap::complex_basis::cartesian<double> cart_3(pol_3);

        std::cout << "pol_1: " << to_string(pol_1) << "\n" << "cart_1: " << to_string(cart_1) << "\n\n";
        REQUIRE(cart_1.covariance()[0][0] != pol_1.covariance()[0][0]);

        // transversity
        const yap::amplitude_basis::transversity<double> t1_pol(pol_1, pol_2, pol_3);
        const yap::amplitude_basis::transversity<double> t1_cart(cart_1, cart_2, cart_3);

        REQUIRE(t1_pol.covariance()[0][0][0][0] != 0);
        REQUIRE(t1_cart.covariance()[0][0][0][0] != 0);
        REQUIRE(t1_pol.covariance()[0][0][0][0] != t1_cart.covariance()[0][0][0][0]);

        // todo: check that amplitudes & covariances actually transform

        const auto c1_cart = yap::amplitude_basis::canonical<double>(t1_cart);
        const auto c1_pol = yap::amplitude_basis::canonical<double>(t1_pol);
        const auto t2 = yap::amplitude_basis::transversity<double>(c1_cart);

        // compare A  pol_x-cart-t1_cart-c1_cart-cart-pol_x_A
        // with    B  pol-       t1_pol -c1_pol -     pol_x_B
        const yap::complex_basis::cartesian<double> cart_1_A(c1_cart.coordinates()[0], c1_cart.covariance()[0][0]);
        const yap::complex_basis::cartesian<double> cart_2_A(c1_cart.coordinates()[1], c1_cart.covariance()[1][1]);
        const yap::complex_basis::cartesian<double> cart_3_A(c1_cart.coordinates()[2], c1_cart.covariance()[2][2]);

        const yap::complex_basis::polar<double> pol_1_A(cart_1_A);
        const yap::complex_basis::polar<double> pol_2_A(cart_2_A);
        const yap::complex_basis::polar<double> pol_3_A(cart_3_A);

        // real and imag really are abs and arg
        const yap::complex_basis::polar<double> pol_1_B(c1_pol.coordinates()[0], c1_pol.covariance()[0][0]);
        const yap::complex_basis::polar<double> pol_2_B(c1_pol.coordinates()[1], c1_pol.covariance()[1][1]);
        const yap::complex_basis::polar<double> pol_3_B(c1_pol.coordinates()[2], c1_pol.covariance()[2][2]);

        std::cout << "cart_1_A: \n" << to_string(cart_1_A) << "\n\n\n";
        std::cout << "pol_1_A/B: \n" << to_string(pol_1_A) << "\n" << to_string(pol_1_B) << "\n\n";
        std::cout << "pol_2_A/B: \n" << to_string(pol_2_A) << "\n" << to_string(pol_2_B) << "\n\n";
        std::cout << "pol_3_A/B: \n" << to_string(pol_3_A) << "\n" << to_string(pol_3_B) << "\n\n";

        REQUIRE( compare_complex_basis(pol_1_A, pol_1_B) );
        REQUIRE( compare_complex_basis(pol_2_A, pol_2_B) );
        REQUIRE( compare_complex_basis(pol_3_A, pol_3_B) );

        REQUIRE( t1_cart.longitudinal()  == t2.longitudinal() );
        REQUIRE( t1_cart.parallel()      == t2.parallel() );
        REQUIRE( t1_cart.perpendicular() == t2.perpendicular() );
        REQUIRE( compare_covariances(t1_cart.covariance(), t2.covariance()) );

        const auto h1 = yap::amplitude_basis::helicity<double>(t1_cart);
        const auto h2 = yap::amplitude_basis::helicity<double>(c1_cart);

        REQUIRE( h1.zero()       == h2.zero()  );
        REQUIRE( h1.plus()       == h2.plus()  );
        REQUIRE( h1.minus()      == h2.minus() );
        REQUIRE( compare_covariances(h1.covariance(), h2.covariance()) );

        const auto t3 = yap::amplitude_basis::transversity<double>(h1);

        REQUIRE( t3.longitudinal()  == t1_cart.longitudinal()  );
        REQUIRE( t3.parallel()      == t1_cart.parallel()      );
        REQUIRE( t3.perpendicular() == t1_cart.perpendicular() );
        REQUIRE( compare_covariances(t3.covariance(), t1_cart.covariance()) );

        const auto c2 = yap::amplitude_basis::canonical<double>(h1);

        REQUIRE( c2[0]           == c1_cart[0] );
        REQUIRE( c2[1]           == c1_cart[1] );
        REQUIRE( c2[2]           == c1_cart[2] );
        REQUIRE( compare_covariances(c2.covariance(), c1_cart.covariance()) );
    }
    SECTION( "non-const" ) {
        // amplitudes in polar coords
        yap::complex_basis::polar<double> pol_1(2.624, yap::rad( 47.), {1.210, 1.951, 0.101});
        yap::complex_basis::polar<double> pol_2(3.157, yap::rad(120.), {2.176, 1.686, 0.153});
        yap::complex_basis::polar<double> pol_3(4.384, yap::rad(163.), {1.782, 2.064, 0.279});

        yap::complex_basis::cartesian<double> cart_1(pol_1);
        yap::complex_basis::cartesian<double> cart_2(pol_2);
        yap::complex_basis::cartesian<double> cart_3(pol_3);

        std::cout << "pol_1: " << to_string(pol_1) << "\n" << "cart_1: " << to_string(cart_1) << "\n\n";
        REQUIRE(cart_1.covariance()[0][0] != pol_1.covariance()[0][0]);

        // transversity
        yap::amplitude_basis::transversity<double> t1_pol(pol_1, pol_2, pol_3);
        yap::amplitude_basis::transversity<double> t1_cart(cart_1, cart_2, cart_3);

        REQUIRE(t1_pol.covariance()[0][0][0][0] != 0);
        REQUIRE(t1_cart.covariance()[0][0][0][0] != 0);
        REQUIRE(t1_pol.covariance()[0][0][0][0] != t1_cart.covariance()[0][0][0][0]);

        auto c1_cart = yap::amplitude_basis::canonical<double>(t1_cart);
        auto c1_pol = yap::amplitude_basis::canonical<double>(t1_pol);
        auto t2 = yap::amplitude_basis::transversity<double>(c1_cart);

        // compare A  pol_x-cart-t1_cart-c1_cart-cart-pol_x_A
        // with    B  pol-       t1_pol -c1_pol -     pol_x_B
        yap::complex_basis::cartesian<double> cart_1_A(c1_cart.coordinates()[0], c1_cart.covariance()[0][0]);
        yap::complex_basis::cartesian<double> cart_2_A(c1_cart.coordinates()[1], c1_cart.covariance()[1][1]);
        yap::complex_basis::cartesian<double> cart_3_A(c1_cart.coordinates()[2], c1_cart.covariance()[2][2]);

        yap::complex_basis::polar<double> pol_1_A(cart_1_A);
        yap::complex_basis::polar<double> pol_2_A(cart_2_A);
        yap::complex_basis::polar<double> pol_3_A(cart_3_A);

        // real and imag really are abs and arg
        yap::complex_basis::polar<double> pol_1_B(c1_pol.coordinates()[0], c1_pol.covariance()[0][0]);
        yap::complex_basis::polar<double> pol_2_B(c1_pol.coordinates()[1], c1_pol.covariance()[1][1]);
        yap::complex_basis::polar<double> pol_3_B(c1_pol.coordinates()[2], c1_pol.covariance()[2][2]);

        std::cout << "cart_1_A: \n" << to_string(cart_1_A) << "\n\n\n";
        std::cout << "pol_1_A/B: \n" << to_string(pol_1_A) << "\n" << to_string(pol_1_B) << "\n\n";
        std::cout << "pol_2_A/B: \n" << to_string(pol_2_A) << "\n" << to_string(pol_2_B) << "\n\n";
        std::cout << "pol_3_A/B: \n" << to_string(pol_3_A) << "\n" << to_string(pol_3_B) << "\n\n";

        REQUIRE( compare_complex_basis(pol_1_A, pol_1_B) );
        REQUIRE( compare_complex_basis(pol_2_A, pol_2_B) );
        REQUIRE( compare_complex_basis(pol_3_A, pol_3_B) );

        REQUIRE( t1_cart.longitudinal()  == t2.longitudinal() );
        REQUIRE( t1_cart.parallel()      == t2.parallel() );
        REQUIRE( t1_cart.perpendicular() == t2.perpendicular() );
        REQUIRE( compare_covariances(t1_cart.covariance(), t2.covariance()) );

        auto h1 = yap::amplitude_basis::helicity<double>(t1_cart);
        auto h2 = yap::amplitude_basis::helicity<double>(c1_cart);

        REQUIRE( h1.zero()       == h2.zero()  );
        REQUIRE( h1.plus()       == h2.plus()  );
        REQUIRE( h1.minus()      == h2.minus() );
        REQUIRE( compare_covariances(h1.covariance(), h2.covariance()) );

        auto t3 = yap::amplitude_basis::transversity<double>(h1);

        REQUIRE( t3.longitudinal()  == t1_cart.longitudinal()  );
        REQUIRE( t3.parallel()      == t1_cart.parallel()      );
        REQUIRE( t3.perpendicular() == t1_cart.perpendicular() );
        REQUIRE( compare_covariances(t3.covariance(), t1_cart.covariance()) );

        auto c2 = yap::amplitude_basis::canonical<double>(h1);

        REQUIRE( c2[0]           == c1_cart[0] );
        REQUIRE( c2[1]           == c1_cart[1] );
        REQUIRE( c2[2]           == c1_cart[2] );
        REQUIRE( compare_covariances(c2.covariance(), c1_cart.covariance()) );
    }
}
