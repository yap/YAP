#include "RelativisticBreitWigner.h"

#include "BlattWeisskopf.h"
#include "CachedValue.h"
#include "Constants.h"
#include "DataPartition.h"
#include "DecayChannel.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "Parameter.h"
#include "Resonance.h"

namespace yap {

//-------------------------
void RelativisticBreitWigner::addDecayChannel(std::shared_ptr<DecayChannel> c)
{
    auto it = resonance()->blattWeisskopfs().find(resonance()->quantumNumbers().twoJ() / 2);
    if (it == resonance()->blattWeisskopfs().end())
        throw exceptions::Exception("Could not store and use Blatt-Weisskopf barrier factor in Resonace_",
                                    "RelativisticBreitWigner::addDecayChannel");
    BlattWeisskopf_ = it->second;
    // add parameters of BlattWeisskopf to this object
    for (const auto& p : BlattWeisskopf_->parameters())
        addParameter(p);
    addParameter(resonance()->radialSize());
    addParameter(mass());
}

//-------------------------
void RelativisticBreitWigner::checkDecayChannel(const DecayChannel& c) const
{
    // check channel is only to two daughters
    if (c.daughters().size() != 2)
        throw exceptions::Exception("Wrong number of daughters (" + std::to_string(c.daughters().size()) + " != 2",
                                    "RelativisticBreitWigner::checkDecayChannel");

    // check channel is to spin-0 particles
    auto sv = spins(c.daughters());
    if (std::any_of(sv.begin(), sv.end(), std::bind(std::not_equal_to<unsigned>(), std::placeholders::_1, 0)))
        throw exceptions::Exception("Non-spin-0 daughter", "RelativisticBreitWigner::checkDecayChannel");
}

//-------------------------
void RelativisticBreitWigner::calculateT(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    /////////////////////////
    // common factors:

    // radial size squared
    double r2 = pow(resonance()->radialSize()->value(), 2);

    // squared resonance mass
    double m2_R = pow(mass()->value(), 2);

    // i * mass^2 * nominal width
    auto im2w_R = Complex_i * m2_R * width()->value();

    // // nominal mass * nominal width
    // double mG = mass()->value() * width()->value();

    // J + 1/2
    unsigned twoLp1 = BlattWeisskopf_->L() + 1;

    /////////////////////////

    // calculate BlattWeisskopf
    BlattWeisskopf_->calculate(D);

    // T := 1 / (M^2 - m^2 - i * M * Gamma)
    for (auto& d : D) {

        // retrieve masses
        double m_ab = model()->fourMomenta()->m(d, pc);
        double m_a  = model()->fourMomenta()->m(d, pc->daughters()[0]);
        double m_b  = model()->fourMomenta()->m(d, pc->daughters()[1]);

        // calculate nominal breakup momentum
        double q2_nomi = squared_breakup_momentum(m2_R, m_a, m_b);

        // retrieve measured breakup momentum
        double q2_meas = model()->measuredBreakupMomenta()->q2(d, pc);

        double Q = sqrt(q2_meas / q2_nomi);

        auto imw = im2w_R / m_ab * pow(Q, twoLp1) * pow(BlattWeisskopf_->value(d, pc), 2) / squared_barrier_factor(BlattWeisskopf_->L(), q2_nomi * r2);

        T()->setValue(1. / (m2_R - pow(m_ab, 2) - imw), d, si, D);

    }
}

}




