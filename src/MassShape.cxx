#include "MassShape.cxx"

#include "Constants.h"

namespace yap {

//-------------------------
MassShape::MassShape() :
    DataAccessor()
{
}

//-------------------------
MassShape::MassShape(const MassShape& other) :
    DataAccessor(other),
    Parameters_(other.Parameters_)
{
}

//-------------------------
MassShape::MassShape(MassShape&& other) :
    DataAccessor(other),
    Parameters_(std::move(other.Parameters_))
{
}

//-------------------------
MassShape::MassShape(const MassShape& rhs)
{
    MassShape temp(rhs);
    std::swap(*this, temp);
    return *this;
}

//-------------------------
MassShape& MassShape::operator=(MassShape&& rhs)
{
    DataAccessor::operator=(std::move(rhs));
    Parameters_ = std::move(rhs.Parameters_);
    return *this;
}

//-------------------------
Amp MassShape::amplitude(DataPoint& d)
{
    return Complex_0;
}

}
