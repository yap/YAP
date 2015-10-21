#include "CachedDataValue.h"

#include "DataAccessor.h"
#include "logging.h"

namespace yap {

//-------------------------
CachedDataValue::CachedDataValue(DataAccessor* owner, unsigned size,
                                 std::vector<std::shared_ptr<ComplexParameter> > ParametersItDependsOn,
                                 std::vector<std::shared_ptr<CachedDataValue> > CachedDataValuesItDependsOn) :
    Owner_(owner),
    Position_(-1),
    Size_(size),
    CalculationStatus_(1, std::vector<CalculationStatus>()),
    VariableStatus_(1, std::vector<VariableStatus>())
{
    if (Size_ == 0)
        LOG(ERROR) << "CachedDataValue::CachedDataValue - Size_ is zero!";
    else if (Owner_) {
        Owner_->addCachedDataValue(this);
        // set position to end of owner's current storage
        Position_ = Owner_->size();
        // increase owner's storage to accomodate cached value
        Owner_->increaseSize(Size_);
    }

    addDependencies(ParametersItDependsOn);
    addDependencies(CachedDataValuesItDependsOn);
}

//-------------------------
CalculationStatus CachedDataValue::calculationStatus(std::shared_ptr<const ParticleCombination> pc, unsigned symmetrizationIndex, unsigned dataPartitionIndex)
{
    // if uncalculated, return without further checking
#ifdef ELPP_DISABLE_DEBUG_LOGS
    if (CalculationStatus_[dataPartitionIndex][symmetrizationIndex] == kUncalculated)
        return kUncalculated;
#else
    if (CalculationStatus_.at(dataPartitionIndex).at(symmetrizationIndex) == kUncalculated)
        return kUncalculated;
#endif

    // else check if any dependencies are changed
    for (auto& c : CachedDataValuesItDependsOn_) {
        // get symmetrization index for particular cached data value,
        // checking if owner is different and grabbing if needed

        LOG(ERROR) << "owner " << c->owner();
        if (! c->owner()->hasSymmetrizationIndex(pc)) {
            LOG(ERROR) << "Error, owner  does not have " << std::string(*pc);
            LOG(ERROR) << "owner " << typeid(*c->owner()).name();
        }

        unsigned symInd = (c->owner() == Owner_) ? symmetrizationIndex : c->owner()->symmetrizationIndex(pc);
        if (c->variableStatus(symInd, dataPartitionIndex) == kChanged) {
            // if dependency has changed, update to uncalculated and return
            setCalculationStatus(kUncalculated, symmetrizationIndex, dataPartitionIndex);
            return kUncalculated;
        }
    }
    for (auto& p : ParametersItDependsOn_) {
        if (p->variableStatus() == kChanged) {
            // if so, update to uncalculated and return
            setCalculationStatus(kUncalculated, symmetrizationIndex, dataPartitionIndex);
            return kUncalculated;
        }
    }

    // else return calculated
    return kCalculated;
}

//-------------------------
void CachedDataValue::setNumberOfSymmetrizations(unsigned n)
{
    for (auto& v : CalculationStatus_) {
        v.resize(n, kUncalculated);
    }
    for (auto& v : VariableStatus_) {
        v.resize(n, kChanged);
    }
}

//-------------------------
void CachedDataValue::setNumberOfDataPartitions(unsigned n)
{
    if (n == 0) {
        LOG(ERROR) << "CachedDataValue::setNumberOfDataPartitions - n cannot be zero.";
        return;
    }
    CalculationStatus_.resize(n, CalculationStatus_.front());
    VariableStatus_.resize(n, VariableStatus_.front());
}

//-------------------------
void RealCachedDataValue::setValue(double val, DataPoint& d, unsigned symmetrizationIndex, unsigned dataPartitionIndex)
{
    if (val != CachedDataValue::value(0, d, symmetrizationIndex)) {
        CachedDataValue::setValue(0, val, d, symmetrizationIndex);
        setVariableStatus(kChanged, symmetrizationIndex, dataPartitionIndex);
    }
    setCalculationStatus(kCalculated, symmetrizationIndex, dataPartitionIndex);
}

//-------------------------
void ComplexCachedDataValue::setValue(double val_re, double val_im, DataPoint& d, unsigned symmetrizationIndex, unsigned dataPartitionIndex)
{
    if (val_re != CachedDataValue::value(0, d, symmetrizationIndex)) {
        CachedDataValue::setValue(0, val_re, d, symmetrizationIndex);
        setVariableStatus(kChanged, symmetrizationIndex, dataPartitionIndex);
    }

    if (val_im != CachedDataValue::value(1, d, symmetrizationIndex)) {
        CachedDataValue::setValue(1, val_im, d, symmetrizationIndex);
        setVariableStatus(kChanged, symmetrizationIndex, dataPartitionIndex);
    }

    setCalculationStatus(kCalculated, symmetrizationIndex, dataPartitionIndex);
}


}
