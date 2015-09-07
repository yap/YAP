#include "BreitWigner.h"

#include "Constants.h"
#include "logging.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double mass, double width) :
    MassShape()
{
    ParameterSet::operator=({mass, width});
}

//-------------------------
Amp BreitWigner::amplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    unsigned sym_index = SymmetrizationIndices_[pc];

    std::vector<double> data = d[index()][sym_index];

    if (CalculationStatuses_[sym_index] == kUncalculated) {
        if (CalcStatus_ == kUncalculated)
            M2iMG_ = Amp(mass() * mass(), -mass() * width());
        /// \todo Create InvariantMassSquared manager inheriting from DataAccessor
        /// for storing m^2 for all relevant particle combinations
        Amp a = Complex_1;
        // Amp a = 1. / (M2iMG_ - Amp(d.invariantMassSquared(pc), 0));
        data[0] = real(a);
        data[1] = imag(a);
        return a;
    }

    return Amp(data[0], data[1]);
}

//-------------------------
Amp BreitWigner::amplitude(double s)
{
    if (CalcStatus_ == kUncalculated)
        M2iMG_ = Amp(mass() * mass(), -mass() * width());
    return 1. / (M2iMG_ - Amp(s, 0));
}

//-------------------------
bool BreitWigner::consistent() const
{
    bool consistent = true;
    if (mass() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - mass <= 0";
        consistent = false;
    }
    if (width() <= 0) {
        LOG(ERROR) << "BreitWigner::consistent() - width <= 0";
        consistent = false;
    }

    return consistent;
}

}




