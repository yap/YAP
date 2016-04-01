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

#include <Constants.h>
#include <CoordinateSystem.h>
#include <HelicityAngles.h>
#include <logging.h>
#include <Rotation.h>

TEST_CASE( "helicityFrame" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    SECTION ( "zero_vector" ) {

        REQUIRE( yap::helicityFrame(yap::ThreeVector_0, yap::ThreeAxes) == yap::ThreeAxes );

    }

    SECTION ( "z" ) {

        REQUIRE( yap::helicityFrame(yap::ThreeAxes[2], yap::ThreeAxes) == yap::ThreeAxes );

    }

    SECTION ( "minus_z" ) {

        auto C = yap::helicityFrame(-yap::ThreeAxes[2], yap::ThreeAxes);

        REQUIRE( C[2] == Catch::Detail::IApprox<yap::ThreeVector<double> >(-yap::ThreeAxes[2]) );
        REQUIRE( C[1] == Catch::Detail::IApprox<yap::ThreeVector<double> >( yap::ThreeAxes[1]) );
        REQUIRE( C[0] == Catch::Detail::IApprox<yap::ThreeVector<double> >(-yap::ThreeAxes[0]) );

    }

    SECTION ( "rotate_Y" ) {

        // rotate by 45 degrees
        auto R = 45 * yap::rad_per_deg<double>();

        // rotate Z around Y by R
        auto P = yap::rotation(yap::ThreeAxes[1], R) * yap::ThreeAxes[2];

        // calculate new helicity frame
        auto C = yap::helicityFrame(P, yap::ThreeAxes);

        // check Z axis of new frame is P
        REQUIRE( C[2] == P );

        // check Y axis is unchanged
        REQUIRE (C[1] == yap::ThreeAxes[1]);

        // check X is rotated by R
        REQUIRE (angle(C[0], yap::ThreeAxes[0]) == Approx(R));

    }

    SECTION ( "rotate_X" ) {

        // rotate by 19.5 degrees
        auto R = 19.5 * yap::rad_per_deg<double>();

        // rotate Z around X by R
        auto P = yap::rotation(yap::ThreeAxes[0], R) * yap::ThreeAxes[2];

        // calculate new helicity frame
        auto C = yap::helicityFrame(P, yap::ThreeAxes);

        // check Z axis of new frame is P
        REQUIRE( C[2] == P );

        // check Y axis is old X axis
        REQUIRE (C[1] == Catch::Detail::IApprox<yap::ThreeVector<double> >(yap::signum(R) * yap::ThreeAxes[0]));
        REQUIRE (C[1] == yap::signum(R) * yap::ThreeAxes[0]);

        // check X is rotated by R from old -Y
        REQUIRE (angle(C[0], -yap::ThreeAxes[1]) == Approx(R));

    }

}
