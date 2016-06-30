#include "CachedDataValue.h"

#include "CalculationStatus.h"
#include "DataAccessor.h"
#include "DataPoint.h"
#include "Exceptions.h"
#include "FourVector.h"
#include "logging.h"
#include "StatusManager.h"
#include "VariableStatus.h"

namespace yap {

//-------------------------
CachedDataValue::Status::Status() :
    Calculation(CalculationStatus::uncalculated),
    Variable(VariableStatus::changed)
{}

CachedDataValue::Status& CachedDataValue::Status::operator=(const VariableStatus& s)
{
    if (Variable != VariableStatus::fixed)
        Variable = s;
    return *this;
}

//-------------------------
std::string to_string(const CachedDataValue::Status& S)
{
    return to_string(S.Calculation) + ", " + to_string(S.Variable);
}

//-------------------------
CachedDataValue::CachedDataValue(unsigned size, ParameterSet pars, CachedDataValueSet vals) :
    std::enable_shared_from_this<CachedDataValue>(),
    Owner_(nullptr),
    Index_(-1),
    Position_(-1),
    Size_(size)
{
    if (Size_ == 0)
        throw exceptions::Exception("zero size", "CachedDataValue::CachedDataValue");
}

//-------------------------
void CachedDataValue::setValue(unsigned index, double val, DataPoint& d, unsigned sym_index) const
{
#ifdef ELPP_DISABLE_DEBUG_LOGS
    d.Data_[Owner_->index()][sym_index][Position_ + index] = val;
#else
    d.Data_.at(Owner_->index()).at(sym_index).at(Position_ + index) = val;
#endif
}

//-------------------------
void CachedDataValue::setDataAccessor(DataAccessor* owner)
{
    if (Owner_ != nullptr)
        throw exceptions::Exception("Owner_ has already been set", "CachedDataValue");

    Owner_ = owner;

    if (Owner_ == nullptr)
        throw exceptions::Exception("Owner_ is nullptr", "CachedDataValue");

    // add self to owner
    Owner_->addCachedDataValue(shared_from_this());

}

//-------------------------
std::shared_ptr<RealCachedDataValue> RealCachedDataValue::create(DataAccessor* da, ParameterSet pars, CachedDataValueSet vals)
{
    auto c = std::shared_ptr<RealCachedDataValue>(new RealCachedDataValue(pars, vals));
    c->setDataAccessor(da);
    return c;
}

//-------------------------
void RealCachedDataValue::setValue(double val, DataPoint& d, unsigned sym_index, StatusManager& sm) const
{
    if (val != CachedDataValue::value(0, d, sym_index)) {
        CachedDataValue::setValue(0, val, d, sym_index);
        sm.status(owner()->index(), index(), sym_index) = VariableStatus::changed;
    }

    sm.status(owner()->index(), index(), sym_index) = CalculationStatus::calculated;
}

//-------------------------
std::shared_ptr<ComplexCachedDataValue> ComplexCachedDataValue::create(DataAccessor* da, ParameterSet pars, CachedDataValueSet vals)
{
    auto c = std::shared_ptr<ComplexCachedDataValue>(new ComplexCachedDataValue(pars, vals));
    c->setDataAccessor(da);
    return c;
}

//-------------------------
void ComplexCachedDataValue::setValue(double val_re, double val_im, DataPoint& d, unsigned sym_index, StatusManager& sm) const
{
    if (val_re != CachedDataValue::value(0, d, sym_index)) {
        CachedDataValue::setValue(0, val_re, d, sym_index);
        sm.status(owner()->index(), index(), sym_index) = VariableStatus::changed;
    }

    if (val_im != CachedDataValue::value(1, d, sym_index)) {
        CachedDataValue::setValue(1, val_im, d, sym_index);
        sm.status(owner()->index(), index(), sym_index) = VariableStatus::changed;
    }

    sm.status(owner()->index(), index(), sym_index) = CalculationStatus::calculated;
}

//-------------------------
std::shared_ptr<FourVectorCachedDataValue> FourVectorCachedDataValue::create(DataAccessor* da, ParameterSet pars, CachedDataValueSet vals)
{
    auto c = std::shared_ptr<FourVectorCachedDataValue>(new FourVectorCachedDataValue(pars, vals));
    c->setDataAccessor(da);
    return c;
}

//-------------------------
const FourVector<double> FourVectorCachedDataValue::value(const DataPoint& d, unsigned  sym_index) const
{
    return FourVector<double>( {CachedDataValue::value(0, d, sym_index),
                                CachedDataValue::value(1, d, sym_index),
                                CachedDataValue::value(2, d, sym_index),
                                CachedDataValue::value(3, d, sym_index)
                               });
}

//-------------------------
void FourVectorCachedDataValue::setValue(const FourVector<double>& val, DataPoint& d, unsigned sym_index) const
{
    for (size_t i = 0; i < val.size(); ++i)
        CachedDataValue::setValue(i, val[i], d, sym_index);
}

//-------------------------
void FourVectorCachedDataValue::setValue(const FourVector<double>& val, DataPoint& d, unsigned sym_index, StatusManager& sm) const
{
    for (size_t i = 0; i < val.size(); ++i) {
        if (val[i] != CachedDataValue::value(i, d, sym_index)) {
            CachedDataValue::setValue(i, val[i], d, sym_index);
            sm.status(owner()->index(), index(), sym_index) = VariableStatus::changed;
        }
    }
    sm.status(owner()->index(), index(), sym_index) = CalculationStatus::calculated;
}

}
