/**
 *
 *    @file  test_ParameterOperations.cxx
 *   @brief  
 *
 *    @date  11/29/16
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include <catch.hpp>
#include <catch_capprox.hpp>

#include <BreitWigner.h>
#include <DataSet.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleTable.h>
#include <ZemachFormalism.h>

#include <cmath>
#include <complex>

TEST_CASE( "ParameterOperations" )
{
    // -- Check operations for real parameters ----------- //
    yap::RealParameter rp(1);

    rp *= 2;
    REQUIRE( rp.value() == 2 );

    rp /= 2;
    REQUIRE( rp.value() == 1 );

    rp += .5;
    REQUIRE( rp.value() == 1.5 );

    rp -= .5;
    REQUIRE( rp.value() == 1.0 );

    // -- Check operations for complex parameters -------- //
    yap::ComplexParameter cp(std::polar<double>(2, yap::rad(30.)));

    cp *= std::polar<double>(2, yap::rad(15.));
    REQUIRE( cp.value() == std::polar<double>(4, yap::rad(45.)) );

    cp /= std::polar<double>(4 / std::sqrt(2), yap::rad(0.));
    std::cout << "Real part of cp = " << std::setprecision(20)
                                      << real(cp.value())      << std::endl;
    std::cout << "Imag part of cp = " << std::setprecision(20)
                                      << imag(cp.value())      << std::endl;
    
    REQUIRE( imag(cp.value()) == Approx(1.) );
    REQUIRE( real(cp.value()) == Approx(1.) );

    // reset cp to (1,1)
    cp.setValue(std::complex<double>(1, 1));

    cp += std::complex<double>(.5, 1);
    REQUIRE( cp.value() == std::complex<double>(1.5, 2) );

    cp -= std::complex<double>(0, .5);
    REQUIRE( cp.value() == std::complex<double>(1.5, 1.5) );
}


