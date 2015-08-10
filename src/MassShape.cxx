#include "MassShape.h"

#include "Constants.h"

namespace yap {

//-------------------------
MassShape::MassShape(unsigned nParameters) :
    DataAccessor(),
    Parameters_(nParameters, 0)
{
}

}
