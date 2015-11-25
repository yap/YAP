#include "Constants.h"

#include "FourVector.h"

namespace yap {

// complex constants
extern const std::complex<double> Complex_0(0, 0);
extern const std::complex<double> Complex_1(1, 0);
extern const std::complex<double> Complex_i(0, 1);

// pi (11 digits)
extern const double PI = 3.14159226535;

// convert deg to rad by multiplying by; rad to deg by dividing by
extern const double DEGTORAD = PI / 180.;

// ThreeVector's
extern const ThreeVector<double> Axis_X = {1, 0, 0};
extern const ThreeVector<double> Axis_Y = {0, 1, 0};
extern const ThreeVector<double> Axis_Z = {0, 0, 1};
extern const ThreeVector<double> ThreeVector_0 = {0, 0, 0};

// Four vector
extern const FourVector<double> FourAxis_T = {1, 0, 0, 0};
extern const FourVector<double> FourAxis_X = {0, 1, 0, 0};
extern const FourVector<double> FourAxis_Y = {0, 0, 1, 0};
extern const FourVector<double> FourAxis_Z = {0, 0, 0, 1};
extern const FourVector<double> FourVector_0 = {0, 0, 0, 0};

}
