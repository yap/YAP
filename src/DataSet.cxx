#include "DataSet.h"

#include "DataPoint.h"

namespace yap {

//-------------------------
bool DataSet::consistent(const DataPoint& d) const
{
    return empty() or equalStructure(front(), d);
}

}
