#include "Constants.h"

namespace yap {

/// complex zero
extern const std::complex<double> Complex_0(0, 0);

/// complex one
extern const std::complex<double> Complex_1(1, 0);

/// complex i
extern const std::complex<double> Complex_i(0, 1);

/// pi (11 digits)
extern const double PI = 3.14159226535;

/// convert deg to rad by multiplying by; rad to deg by dividing by
extern const double DEGTORAD = PI / 180.;

}
