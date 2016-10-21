#include "CachedValue.h"

#include "CalculationStatus.h"
#include "DataAccessor.h"
#include "DataPoint.h"
#include "Exceptions.h"
#include "FourVector.h"
#include "StatusManager.h"
#include "VariableStatus.h"

namespace yap {

//-------------------------
CachedValue::Status::Status() :
    Calculation(CalculationStatus::uncalculated),
    Variable(VariableStatus::changed)
{}

CachedValue::Status& CachedValue::Status::operator=(const VariableStatus& s)
{
    if (Variable != VariableStatus::fixed)
        Variable = s;
    return *this;
}

//-------------------------
std::string to_string(const CachedValue::Status& S)
{
    return to_string(S.Calculation) + ", " + to_string(S.Variable);
}

//-------------------------
CachedValue::CachedValue(unsigned size, DataAccessor& da) :
    std::enable_shared_from_this<CachedValue>(),
    Owner_(&da),
    Index_(-1),
    Position_(-1),
    Size_(size)
{
    if (Size_ == 0)
        throw exceptions::Exception("zero size", "CachedValue::CachedValue");
}

//-------------------------
void CachedValue::addToDataAccessor()
{
    // add self to owner
    Owner_->addCachedValue(*this);
}

//-------------------------
std::shared_ptr<RealCachedValue> RealCachedValue::create(DataAccessor& da)
{
    auto c = std::shared_ptr<RealCachedValue>(new RealCachedValue(da));
    c->addToDataAccessor();
    return c;
}

//-------------------------
void RealCachedValue::setValue(double val, DataPoint& d, unsigned sym_index, StatusManager& sm) const
{
    if (val != CachedValue::value(0, d, sym_index)) {
        CachedValue::setValue(0, val, d, sym_index);
        sm.status(owner()->index(), index(), sym_index) = VariableStatus::changed;
    }

    sm.status(owner()->index(), index(), sym_index) = CalculationStatus::calculated;
}

//-------------------------
std::shared_ptr<ComplexCachedValue> ComplexCachedValue::create(DataAccessor& da)
{
    auto c = std::shared_ptr<ComplexCachedValue>(new ComplexCachedValue(da));
    c->addToDataAccessor();
    return c;
}

//-------------------------
void ComplexCachedValue::setValue(double val_re, double val_im, DataPoint& d, unsigned sym_index, StatusManager& sm) const
{
    if (val_re != CachedValue::value(0, d, sym_index)) {
        CachedValue::setValue(0, val_re, d, sym_index);
        sm.status(owner()->index(), index(), sym_index) = VariableStatus::changed;
    }

    if (val_im != CachedValue::value(1, d, sym_index)) {
        CachedValue::setValue(1, val_im, d, sym_index);
        sm.status(owner()->index(), index(), sym_index) = VariableStatus::changed;
    }

    sm.status(owner()->index(), index(), sym_index) = CalculationStatus::calculated;
}

//-------------------------
std::shared_ptr<FourVectorCachedValue> FourVectorCachedValue::create(DataAccessor& da)
{
    auto c = std::shared_ptr<FourVectorCachedValue>(new FourVectorCachedValue(da));
    c->addToDataAccessor();
    return c;
}

//-------------------------
const FourVector<double> FourVectorCachedValue::value(const DataPoint& d, unsigned  sym_index) const
{
    return FourVector<double>( {CachedValue::value(0, d, sym_index),
                                CachedValue::value(1, d, sym_index),
                                CachedValue::value(2, d, sym_index),
                                CachedValue::value(3, d, sym_index)
                               });
}

//-------------------------
void FourVectorCachedValue::setValue(const FourVector<double>& val, DataPoint& d, unsigned sym_index) const
{
    for (size_t i = 0; i < val.size(); ++i)
        CachedValue::setValue(i, val[i], d, sym_index);
}

//-------------------------
void FourVectorCachedValue::setValue(const FourVector<double>& val, DataPoint& d, unsigned sym_index, StatusManager& sm) const
{
    for (size_t i = 0; i < val.size(); ++i) {
        if (val[i] != CachedValue::value(i, d, sym_index)) {
            CachedValue::setValue(i, val[i], d, sym_index);
            sm.status(owner()->index(), index(), sym_index) = VariableStatus::changed;
        }
    }
    sm.status(owner()->index(), index(), sym_index) = CalculationStatus::calculated;
}

}
