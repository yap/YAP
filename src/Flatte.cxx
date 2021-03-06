#include "Flatte.h"

#include "CalculationStatus.h"
#include "DataPoint.h"
#include "DataPartition.h"
#include "DecayChannel.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleTable.h"

namespace yap {

//-------------------------
FlatteChannel::FlatteChannel(std::shared_ptr<NonnegativeRealParameter> coupling, FinalStateParticle& A, FinalStateParticle& B) :
    Coupling(coupling),
    Particles({std::static_pointer_cast<FinalStateParticle>(A.shared_from_this()),
                std::static_pointer_cast<FinalStateParticle>(B.shared_from_this())})
{
}

//-------------------------
FlatteChannel::FlatteChannel(double coupling, FinalStateParticle& A, FinalStateParticle& B) :
    FlatteChannel(std::make_shared<NonnegativeRealParameter>(coupling), A, B)
{
}

//-------------------------
FlatteChannel::FlatteChannel(double coupling, const ParticleTableEntry& a, const ParticleTableEntry& b) :
    FlatteChannel(coupling, *FinalStateParticle::create(a), *FinalStateParticle::create(b))
{
}

//-------------------------
Flatte::Flatte(double m) :
    MassShapeWithNominalMass(m),
    T_(ComplexCachedValue::create(*this))
{
}

//-------------------------
Flatte::Flatte(const ParticleTableEntry& pde) :
    MassShapeWithNominalMass(pde),
    T_(ComplexCachedValue::create(*this))
{
}

//-------------------------
void Flatte::add(FlatteChannel fc)
{
    if (!fc.Coupling)
        throw exceptions::Exception("Coupling is unset", "Flatte::add");
    if (!fc.Particles[0] or !fc.Particles[1])
        throw exceptions::Exception("Particles unset", "Flatte::add");
    // check if channel already contained
    for (const auto& FC : FlatteChannels_)
        if (FC.Particles == fc.Particles or (FC.Particles[0] == fc.Particles[1] and FC.Particles[1] == fc.Particles[0]))
            throw exceptions::Exception("Channel already held", "Flatte::add");

    FlatteChannels_.push_back(fc);

    addParameter(FlatteChannels_.back().Coupling);
}

//-------------------------
void Flatte::checkDecayChannel(const DecayChannel& c) const
{
    // check Flatte has channels
    if (FlatteChannels_.empty())
        throw exceptions::Exception("Add FlatteChannels to Flatte before adding DecayChannel to DecayingParticle", "Flatte::checkDecayChannel");

    // check channel has right number of daughters
    if (c.daughters().size() != FlatteChannels_[0].Particles.size())
        throw exceptions::Exception("Wrong number of daughters (" + std::to_string(c.daughters().size())
                                    + " != " + std::to_string(FlatteChannels_[0].Particles.size()),
                                    "Flatte::checkDecayChannel");

    // check that Flatte has FlatteChannel with correct particle content
    for (const auto& fc : FlatteChannels_)
        if ((fc.Particles[0] == c.daughters()[0] and fc.Particles[1] == c.daughters()[1])
            or
            (fc.Particles[0] == c.daughters()[1] and fc.Particles[1] == c.daughters()[0]))
            // if found, return before exception can be thrown
            return;
    throw exceptions::Exception("Flatte doesn't contain channel for " + to_string(c), "Flatte::checkDecayChannel");
}
    
//-------------------------
void Flatte::updateCalculationStatus(StatusManager& D) const
{
    if (status() == VariableStatus::changed)
        D.set(*T_, CalculationStatus::uncalculated);
}

//-------------------------
const std::complex<double> Flatte::value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    return T_->value(d, symmetrizationIndex(pc));
}
    
//-------------------------
void Flatte::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T_, si) != CalculationStatus::uncalculated)
        return;

    /////////////////////////
    // common factors:

    // mass^2
    const auto m2 = pow(mass()->value(), 2);

    // width at nominal mass
    std::complex<double> w = 0;
    for (const auto& fc : FlatteChannels_)
        w += fc.Coupling->value() * std::sqrt(std::complex<double>(measured_breakup_momenta::q2(m2, fc.Particles[0]->mass(), fc.Particles[1]->mass())));

    // normalization factor
    auto w_o_m = 2. * w / mass()->value();

    /////////////////////////

    for (auto& d : D) {
        
        auto s = model()->fourMomenta()->m2(d, pc);

        // calculate width term := sum of coupling * complex-breakup-momentum
        std::complex<double> ws = 0;
        for (const auto& fc : FlatteChannels_)
            ws += fc.Coupling->value() * std::sqrt(std::complex<double>(measured_breakup_momenta::q2(s, fc.Particles[0]->mass(), fc.Particles[1]->mass())));

        // T = 1 / (M^2 - m^2 - width-term)
        T_->setValue(w_o_m / (m2 - s - 1_i * 2. * ws / sqrt(s)), d, si, D);
    }

    D.status(*T_, si) = CalculationStatus::calculated;
}

}




