#include "Integrator.h"

#include "DecayTreeVectorIntegral.h"
#include "IntegralElement.h"
#include "ModelIntegral.h"

namespace yap {

//-------------------------
IntegralMap& Integrator::integrals(ModelIntegral& I)
{ return I.Integrals_; }

//-------------------------
RealIntegralElementVector& Integrator::diagonals(DecayTreeVectorIntegral& I)
{ return I.Diagonals_; }

//-------------------------
ComplexIntegralElementMatrix& Integrator::offDiagonals(DecayTreeVectorIntegral& I)
{ return I.OffDiagonals_; }

//-------------------------
DecayTreeVectorIntegral& Integrator::reset(DecayTreeVectorIntegral& I)
{ return I.reset(); }

}
