#include "Filters.h"

#include "container_utils.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "DecayTree.h"
#include "FinalStateParticle.h"
#include "FreeAmplitude.h"
#include "MassShape.h"
#include "MassShapeWithNominalMass.h"
#include "Particle.h"
#include "Resonance.h"
#include "SpinAmplitude.h"

namespace yap {

//-------------------------
const bool filter_free_amplitude::operator()(const DecayTree& dt) const
{
    return by_ptr(*this, dt.freeAmplitude());
}

//-------------------------
const bool filter_parameter::operator()(const FreeAmplitude& fa) const
{
    return operator()(static_cast<const ParameterBase&>(fa));
}

//-------------------------
const bool filter_decay_channel::operator()(const FreeAmplitude& fa) const
{
    return by_ptr(*this, fa.decayChannel());
}

//-------------------------
const bool filter_spin_amplitude::operator()(const FreeAmplitude& fa) const
{
    return by_ptr(*this, fa.spinAmplitude());
}

//-------------------------
const bool filter_decaying_particle::operator()(const Particle& p) const
{
    if (dynamic_cast<const DecayingParticle*>(&p))
        return operator()(static_cast<const DecayingParticle&>(p));
    return false;
}

//-------------------------
const bool filter_resonance::operator()(const DecayingParticle& p) const
{
    if (dynamic_cast<const Resonance*>(&p))
        return operator()(static_cast<const Resonance&>(p));
    return false;
}

//-------------------------
const bool filter_mass_shape::operator()(const Resonance& p) const
{
    if (p.massShape())
        return operator()(*p.massShape());
    return false;
}

//-------------------------
const bool has_decay_channel::operator()(const FreeAmplitude& fa) const
{
    return operator()(fa.decayChannel());
}

//-------------------------
const bool has_decay_channel::operator()(const DecayingParticle& dp) const
{
    return std::any_of(dp.channels().begin(), dp.channels().end(), *this);
}

//-------------------------
const bool has_spin_amplitude::operator()(const FreeAmplitude& fa) const
{
    return operator()(fa.spinAmplitude());
}

//-------------------------
const bool to::operator()(const DecayChannel& dc) const
{
    return orderless_contains(dc.daughters().begin(), dc.daughters().end(), Daughters_.begin(), Daughters_.end());
}

//-------------------------
const bool to::operator()(const DecayingParticle& dp) const
{
    for (const auto& dc : dp.channels())
        if (by_ptr(*this, dc))
            return true;
    return false;
}

//-------------------------
const bool exactly_to::operator()(const DecayChannel& dc) const
{
    return orderless_equal(dc.daughters().begin(), dc.daughters().end(), Daughters_.begin(), Daughters_.end());
}

//-------------------------
const bool l_equals::operator()(const SpinAmplitude& sa) const
{
    return sa.L() == L_;
}

//-------------------------
const bool m_equals::operator()(const FreeAmplitude& fa) const
{
    return fa.twoM() == TwoM_;
}

//-------------------------
const bool has_decay_tree::operator()(const DecayingParticle& dp) const
{
    for (const auto& m_dtv : dp.decayTrees())
        if (std::any_of(m_dtv.second.begin(), m_dtv.second.end(), *this))
            return true;
    return false;
}

//-------------------------
const bool has_free_amplitude::operator()(const DecayTree& dp) const
{
    return operator()(dp.freeAmplitude());
}

//-------------------------
const bool has_free_amplitude::operator()(const DecayingParticle& dp) const
{
    for (const auto& m_dtv : dp.decayTrees())
        for (const auto& dt : m_dtv.second)
            if (by_ptr(*this, dt))
                return true;
    return false;
}

//-------------------------
const bool from::operator()(const Particle& p) const
{
    return to(std::const_pointer_cast<Particle>(p.shared_from_this()))(*DecayingParticle_);
}

//-------------------------
const bool has_mass::operator()(const Particle& p) const
{
    if (dynamic_cast<const DecayingParticle*>(&p))
        return operator()(*static_cast<const DecayingParticle*>(&p));
    if (dynamic_cast<const FinalStateParticle*>(&p))
        return operator()(*static_cast<const FinalStateParticle*>(&p));
    return false;
}

//-------------------------
const bool has_mass::operator()(const MassShape& m) const
{
    if (dynamic_cast<const MassShapeWithNominalMass*>(&m))
        return operator()(*static_cast<const MassShapeWithNominalMass*>(&m));
    return false;
}

//-------------------------
const bool is_fixed::operator()(const ParameterBase& p) const
{
    return p.variableStatus() == VariableStatus::fixed;
}

//-------------------------
const bool is_not_fixed::operator()(const ParameterBase& p) const
{
    return p.variableStatus() != VariableStatus::fixed;
}

}
