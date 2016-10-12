/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/// \file
/// \brief Tests calculation of helicity reference frame

#include <catch.hpp>
#include <catch_iapprox.hpp>

#include <CoordinateSystem.h>
#include <HelicityAngles.h>
#include <logging.h>
#include <MathUtilities.h>
#include <Rotation.h>

TEST_CASE( "helicityFrame" )
{

    yap::CoordinateSystem<double, 3> three_axes({yap::ThreeVector<double>({1., 0., 0.}),
                                                 yap::ThreeVector<double>({0., 1., 0.}),
                                                 yap::ThreeVector<double>({0., 0., 1.})});

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    SECTION ( "zero_vector" ) {

        REQUIRE( yap::helicityFrame(yap::ThreeVector<double>({0, 0, 0}), three_axes) == three_axes );

    }

    SECTION ( "z" ) {

        REQUIRE( yap::helicityFrame(three_axes[2], three_axes) == three_axes );

    }

    SECTION ( "minus_z" ) {

        auto C = yap::helicityFrame(-three_axes[2], three_axes);

        REQUIRE( C[2] == Catch::Detail::IApprox<yap::ThreeVector<double> >(-three_axes[2]) );
        REQUIRE( C[1] == Catch::Detail::IApprox<yap::ThreeVector<double> >( three_axes[1]) );
        REQUIRE( C[0] == Catch::Detail::IApprox<yap::ThreeVector<double> >(-three_axes[0]) );

    }

    SECTION ( "rotate_Y" ) {

        // rotate by 45 degrees
        auto R = yap::rad(45.);

        // rotate Z around Y by R
        auto P = yap::rotation(three_axes[1], R) * three_axes[2];

        // calculate new helicity frame
        auto C = yap::helicityFrame(P, three_axes);

        // check Z axis of new frame is P
        REQUIRE( C[2] == P );

        // check Y axis is unchanged
        REQUIRE (C[1] == three_axes[1]);

        // check X is rotated by R
        REQUIRE (angle(C[0], three_axes[0]) == Approx(R));

    }

    SECTION ( "rotate_X" ) {

        // rotate by 19.5 degrees
        auto R = yap::rad(19.5);

        // rotate Z around X by R
        auto P = yap::rotation(three_axes[0], R) * three_axes[2];

        // calculate new helicity frame
        auto C = yap::helicityFrame(P, three_axes);

        // check Z axis of new frame is P
        REQUIRE( C[2] == P );

        // check Y axis is old X axis
        REQUIRE (C[1] == Catch::Detail::IApprox<yap::ThreeVector<double> >(yap::signum(R) * three_axes[0]));
        REQUIRE (C[1] == yap::signum(R) * three_axes[0]);

        // check X is rotated by R from old -Y
        REQUIRE (angle(C[0], -three_axes[1]) == Approx(R));

    }

}
