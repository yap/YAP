#include "Integrator.h"

#include "DecayTreeVectorIntegral.h"
#include "IntegralElement.h"
#include "ModelIntegral.h"

namespace yap {

//-------------------------
IntegralMap& Integrator::integrals(ModelIntegral& I)
{ return I.Integrals_; }

//-------------------------
DiagonalIntegralMap& Integrator::diagonals(DecayTreeVectorIntegral& I)
{ return I.Diagonals_; }

//-------------------------
OffDiagonalIntegralMap& Integrator::offDiagonals(DecayTreeVectorIntegral& I)
{ return I.OffDiagonals_; }

//-------------------------
DiagonalIntegralMap::mapped_type& Integrator::diagonalComponent(DecayTreeVectorIntegral& I, const DecayTreeVector::value_type& dt)
{ return I.Diagonals_.at(dt); }

//-------------------------
OffDiagonalIntegralMap::mapped_type& Integrator::offDiagonalComponent(DecayTreeVectorIntegral& I, const DecayTreeVector::value_type& i, const DecayTreeVector::value_type& j)
{ return I.OffDiagonals_.at(OffDiagonalIntegralMap::key_type({i, j})); }

}
