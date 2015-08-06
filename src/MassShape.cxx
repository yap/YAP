#include "MassShape.cxx"

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
MassShape& MassShape::operator=(const MassShape& rhs)
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
