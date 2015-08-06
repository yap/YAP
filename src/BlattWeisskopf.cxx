#include "BlattWeisskopf.h"

namespace yap {

//-------------------------
Amp BlattWeisskopf::amplitude(DataPoint& d)
{
    return Amp(1);
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    return true;
}

}

