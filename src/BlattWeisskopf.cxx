#include "BlattWeisskopf.h"

namespace yap {

//-------------------------
BlattWeisskopf::BlattWeisskopf(DecayChannel* decayChannel) :
    AmplitudeComponent(),
    DecayChannel_(decayChannel)
{}

//-------------------------
Amp BlattWeisskopf::amplitude(DataPoint& d)
{
    /// \todo implement
    return Amp(1);
}

//-------------------------
bool BlattWeisskopf::consistent() const
{
    /// \todo implement
    return true;
}

}

