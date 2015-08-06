#include "MassShape.cxx"

//-------------------------
MassShape::MassShape()
{
}

//-------------------------
MassShape::MassShape(const MassShape& other)
    : Parameters_(other.Parameters_)
{
}

//-------------------------
MassShape::MassShape(MassShape&& other)
    : Parameters_(std::move(other.Parameters_))
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
    Parameters_ = std::move(rhs.Parameters_);
    return *this;
}
