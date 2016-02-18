#include "CachedDataValue.h"

#include "DataPoint.h"
#include "Exceptions.h"
#include "logging.h"

namespace yap {

//-------------------------
CachedDataValue::CachedDataValue(unsigned size, ParameterSet pars, CachedDataValueSet vals) :
    std::enable_shared_from_this<CachedDataValue>(),
    Owner_(nullptr),
    Position_(-1),
    Size_(size),
    ParametersItDependsOn_(pars),
    CachedDataValuesItDependsOn_(vals),
    CalculationStatus_(1, std::vector<CalculationStatus>()),
    VariableStatus_(1, std::vector<VariableStatus>())
{
    if (Size_ == 0)
        throw exceptions::Exception("zero size", "CachedDataValue::CachedDataValue");
}

//-------------------------
bool CachedDataValue::consistent() const
{
    return false;
}

//-------------------------
CalculationStatus CachedDataValue::globalCalculationStatus(const std::shared_ptr<ParticleCombination>& pc)
{
    return GlobalCalculationStatus_[Owner_->symmetrizationIndex(pc)];
}

//-------------------------
double CachedDataValue::value(unsigned index, const DataPoint& d, unsigned symmetrizationIndex) const
{
#ifdef ELPP_DISABLE_DEBUG_LOGS
    return d.Data_[Owner_->index()][symmetrizationIndex][Position_ + index];
#else
    return d.Data_.at(Owner_->index()).at(symmetrizationIndex).at(Position_ + index);
#endif
}

//-------------------------
void CachedDataValue::setCalculationStatus(CalculationStatus stat, const std::shared_ptr<ParticleCombination>& pc, unsigned dataPartitionIndex)
{
    CalculationStatus_[dataPartitionIndex][Owner_->symmetrizationIndex(pc)] = stat;
}

//-------------------------
void CachedDataValue::setGlobalCalculationStatus(CalculationStatus stat, const std::shared_ptr<ParticleCombination>& pc)
{
    GlobalCalculationStatus_[Owner_->symmetrizationIndex(pc)] = stat;
}

//-------------------------
void CachedDataValue::setValue(unsigned index, double val, DataPoint& d, unsigned symmetrizationIndex) const
{
    d.Data_[Owner_->index()][symmetrizationIndex][Position_ + index] = val;
}

//-------------------------
void CachedDataValue::updateGlobalCalculationStatus(const std::shared_ptr<ParticleCombination>& pc)
{
    updateGlobalCalculationStatus(pc, Owner_->symmetrizationIndex(pc));
}

//-------------------------
void CachedDataValue::removeDependency(std::shared_ptr<ParameterBase> dep)
{
    // look for parameter in set of dependencies
    auto it = ParametersItDependsOn_.find(dep);
    // if found, erase
    if (it != ParametersItDependsOn_.end())
        ParametersItDependsOn_.erase(it);
}

//-------------------------
void CachedDataValue::updateGlobalCalculationStatus(const std::shared_ptr<ParticleCombination>& pc, unsigned symmetrizationIndex)
{
    DEBUG("CachedDataValue::updateGlobalCalculationStatus - symIndex for " << *pc << " = " << symmetrizationIndex);

    // if CachedDataValue is uncalculated for any of the DataPartitions, set to uncalculated
    // to make sure it will be calculated during the first iteration
    for (auto& v : CalculationStatus_) {
        if (v[symmetrizationIndex] == kUncalculated) {
            GlobalCalculationStatus_[symmetrizationIndex] = kUncalculated;
            DEBUG("kUncalculated (CachedDataValue is uncalculated for some of the DataPartitions)");
            return;
        }
    }

    // check if any dependencies are changed
    for (auto& p : ParametersItDependsOn_) {
        if (p->variableStatus() == kChanged) {
            // if so, update to uncalculated and return
            GlobalCalculationStatus_[symmetrizationIndex] = kUncalculated;
            DEBUG("kUncalculated (ParametersItDependsOn_ are changed)");
            return;
        }
    }

    for (auto& c : CachedDataValuesItDependsOn_) {

        // if the owner does not have the symIndex, there is nothing to check
        if (c->owner() != Owner_ and not c->owner()->hasParticleCombination(pc))
            continue;


        DEBUG(" updateGlobalCalculationStatus of CachedDataValueItDependsOn");
        c->updateGlobalCalculationStatus(pc);
        DEBUG(" done updateGlobalCalculationStatus of CachedDataValueItDependsOn");

        if (c->globalCalculationStatus(pc) == kUncalculated) {
            GlobalCalculationStatus_[symmetrizationIndex] = kUncalculated;
            DEBUG("kUncalculated (globalCalculationStatus is kUncalculated)");
            return;
        }
    }

    for (auto& c : DaughterCachedDataValuesItDependsOn_) {

        const std::shared_ptr<ParticleCombination>& cPc = pc->daughters()[c.second];

        // if the owner does not have the symIndex, there is nothing to check
        if (c.first->owner() != Owner_ and not c.first->owner()->hasParticleCombination(cPc))
            continue;


        DEBUG(" updateGlobalCalculationStatus of DaughterCachedDataValueItDependsOn");
        c.first->updateGlobalCalculationStatus(cPc);
        DEBUG(" done updateGlobalCalculationStatus of DaughterCachedDataValueItDependsOn");

        if (c.first->globalCalculationStatus(cPc) == kUncalculated) {
            GlobalCalculationStatus_[symmetrizationIndex] = kUncalculated;
            DEBUG("kUncalculated (globalCalculationStatus is kUncalculated)");
            return;
        }
    }

    // otherwise nothing has changed and it is calculated
    GlobalCalculationStatus_[symmetrizationIndex] = kCalculated;
    DEBUG("kCalculated (nothing has changed)");
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
    setNumberOfSymmetrizations(Owner_->maxSymmetrizationIndex() + 1);
    Owner_->addCachedDataValue(shared_from_this());

    // set position to end of owner's current storage
    Position_ = Owner_->size();

    // increase owner's storage to accommodate cached value
    Owner_->increaseSize(Size_);
}

//-------------------------
void CachedDataValue::setNumberOfSymmetrizations(unsigned n)
{
    for (auto& v : CalculationStatus_) {
        v.resize(n, kUncalculated);
    }

    GlobalCalculationStatus_.resize(n, kUncalculated);

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
std::shared_ptr<RealCachedDataValue> RealCachedDataValue::create(DataAccessor* da, ParameterSet pars, CachedDataValueSet vals)
{
    auto c = std::shared_ptr<RealCachedDataValue>(new RealCachedDataValue(pars, vals));
    c->setDataAccessor(da);
    return c;
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
std::shared_ptr<ComplexCachedDataValue> ComplexCachedDataValue::create(DataAccessor* da, ParameterSet pars, CachedDataValueSet vals)
{
    auto c = std::shared_ptr<ComplexCachedDataValue>(new ComplexCachedDataValue(pars, vals));
    c->setDataAccessor(da);
    return c;
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

//-------------------------
std::shared_ptr<FourVectorCachedDataValue> FourVectorCachedDataValue::create(DataAccessor* da, ParameterSet pars, CachedDataValueSet vals)
{
    auto c = std::shared_ptr<FourVectorCachedDataValue>(new FourVectorCachedDataValue(pars, vals));
    c->setDataAccessor(da);
    return c;
}

//-------------------------
void FourVectorCachedDataValue::setValue(FourVector<double> val, DataPoint& d, unsigned symmetrizationIndex, unsigned dataPartitionIndex)
{
    for (size_t i = 0; i < val.size(); ++i) {
        if (val[i] != CachedDataValue::value(i, d, symmetrizationIndex)) {
            CachedDataValue::setValue(i, val[i], d, symmetrizationIndex);
            setVariableStatus(kChanged, symmetrizationIndex, dataPartitionIndex);
        }
    }
    setCalculationStatus(kCalculated, symmetrizationIndex, dataPartitionIndex);
}

}
