#include "DataSet.h"

#include "DataPoint.h"
#include "FourMomenta.h"
#include "logging.h"

namespace yap {

//-------------------------
bool DataSet::consistent(const DataPoint& d) const
{
    return empty() or equalStructure(front(), d);
}

}
