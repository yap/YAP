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

#include <Constants.h>
#include <CoordinateSystem.h>
#include <HelicityAngles.h>
#include <logging.h>
#include <Rotation.h>

TEST_CASE( "helicityFrame" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);

    SECTION ( "zero_vector" ) {

        REQUIRE( yap::helicityFrame(yap::ThreeVector_0, yap::ThreeAxes) == yap::ThreeAxes );

    }

    SECTION ( "equal_z" ) {

        REQUIRE( yap::helicityFrame(yap::ThreeAxes[2], yap::ThreeAxes) == yap::ThreeAxes );

    }

    SECTION ( "forty-five_degrees" ) {

        // rotate Z around Y by 45 degrees
        auto P = yap::rotation(yap::ThreeAxes[1], yap::pi<double>() / 4) * yap::ThreeAxes[2];

        // calculate new helicity frame
        auto C = yap::helicityFrame(P, yap::ThreeAxes);

        // check Z axis of new frame is P
        REQUIRE( C[2] == P );

        // check Y axis is unchanged
        REQUIRE (C[1] == yap::ThreeAxes[1]);

        // check X is rotated by 45 degrees
        REQUIRE (angle(C[0], yap::ThreeAxes[0]) == Approx(yap::pi<double>() / 4.));

    }

}
