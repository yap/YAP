#include "DataAccessor.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor() :
    recalculate_(true),
    index_(0)
{
  // assign a running index to this DataAccessor
  static unsigned int index(0);
  index_ = index++;
}

}
