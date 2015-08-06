#include "BreitWigner.h"

#include "Constants.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double mass, double width) :
    MassShape()
{
    setMass(mass);
    setWidth(width);
}

//-------------------------
BreitWigner::BreitWigner(const BreitWigner& other) :
    MassShape(other)
{
}

//-------------------------
BreitWigner::BreitWigner(BreitWigner&& other) :
    MassShape(other)
{
}

//-------------------------
BreitWigner::~BreitWigner()
{
}

//-------------------------
BreitWigner& BreitWigner::operator=(const BreitWigner& rhs)
{
    MassShape::operator=(rhs);
    return *this;
}

//-------------------------
BreitWigner& BreitWigner::operator=(BreitWigner&& rhs)
{
    MassShape::operator=(rhs);
    return *this;
}

//-------------------------
Amp BreitWigner::amplitude(DataPoint& d)
{
    return Complex_0;
}

//-------------------------
Amp BreitWigner::amplitude(double s)
{
    return 1. / Amp(squaredmass(s) - s, mass() * width(s));
}

//-------------------------
bool BreitWigner::consistent() const
{
    return (mass() > 0 && width() > 0);
}

}




