#include "DataAccessor.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor() :
    Recalculate_(true),
    Index_(0)
{
    // assign a running index to this DataAccessor
    // TODO: come up with something smarter
    static unsigned int index(0);
    Index_ = index++;
}

}
