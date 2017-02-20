#include "DalitzPhspIntegral.h"

#include "FastNumericalIntegration/DEIntegrator.h"
#include "FinalStateParticle.h"
#include "logging.h"

#include <assert.h>

namespace yap {

//-------------------------
DalitzIntegrand::DalitzIntegrand(double isp_mass, const FinalStateParticleVector& fsps) :
    M2_(isp_mass*isp_mass)
{
    if (fsps.size() != 3)
        throw exceptions::Exception("Can only calculate Dalitz phasespace volume for 3 final state particles", "DalitzIntegrand::DalitzIntegrand");

    ma2_ = pow(fsps[0]->mass(), 2);
    mb2_ = pow(fsps[1]->mass(), 2);
    mc2_ = pow(fsps[2]->mass(), 2);
}

//-------------------------
double DalitzIntegrand::operator()(double mab2) const
{
    double pb2 = pow(mab2 - ma2_ + mb2_, 2)/(4.*mab2) - mb2_; // p_b^2 = E_b^2 - m_b^2
    double pc2 = pow(M2_ - mab2 - mc2_,  2)/(4.*mab2) - mc2_; // p_c^2 = E_c^2 - m_c^2

    if (pb2 <= 0 or pc2 <= 0)
        return 0;

    return 4. * sqrt(pb2*pc2);
}

//-------------------------
double dalitz_phasespace_volume(double isp_mass, const FinalStateParticleVector& fsps, const double relErr) {
    if (fsps.size() != 3)
        throw exceptions::Exception("Can only calculate Dalitz phasespace volume for 3 final state particles", "dalitz_phasespace_volume");

    const double prefactor = 1. / (isp_mass*isp_mass); // \todo is this correct ???
    //const double prefactor = 1.;

    DalitzIntegrand f(isp_mass, fsps);
    double lowerBound = pow(fsps[0]->mass() + fsps[1]->mass(), 2); // (m_a + m_b)^2
    double upperBound = pow(isp_mass - fsps[2]->mass(), 2); // (M - m_c)^2
    const double volume_upper_bound = prefactor * pow(upperBound - lowerBound, 2);

    double absErr = relErr * 0.5 * volume_upper_bound;

    int evaluations(0);
    double errorEstimate(1.e99); // absolute
    double result(0.);

    for (unsigned i = 0; i < 10; ++i) {
        result = prefactor *DEIntegrator<DalitzIntegrand>::Integrate(f, lowerBound, upperBound, absErr, evaluations, errorEstimate);
        assert (result <= volume_upper_bound);

        DEBUG("Dalitz phasespace volume = " << result << "; calculated with " << evaluations << "evaluations");

        if (errorEstimate/result <= relErr)
            break;

        absErr = 0.8 * std::min(result * relErr, absErr); // decrease absErr
    }

    return result;
}

}
