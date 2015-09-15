#include "DataAccessor.h"

#include "DataPoint.h"
#include "InitialStateParticle.h"
#include "logging.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor(InitialStateParticle* isp, ParticleCombination::Equiv* equiv) :
    InitialStateParticle_(isp),
    Equiv_(equiv),
    Index_(0)
{
    if (isp == nullptr)
        LOG(ERROR) << "DataAccessor: no InitialStateParticle provided!";

    if (this != InitialStateParticle_)
        InitialStateParticle_->addDataAccessor(this);

    // index is set later by InitialStateParticle::setDataAcessorIndices()
    // via InitialStateParticle::prepare()
}

//-------------------------
DataAccessor::DataAccessor(const DataAccessor& other) :
    InitialStateParticle_(other.InitialStateParticle_),
    Equiv_(other.Equiv_),
    SymmetrizationIndices_(other.SymmetrizationIndices_),
    Index_(0)
{
    InitialStateParticle_->addDataAccessor(this);
}

//-------------------------
DataAccessor::~DataAccessor()
{
    if (this != InitialStateParticle_)
        InitialStateParticle_->removeDataAccessor(this);
}

//-------------------------
unsigned DataAccessor::maxSymmetrizationIndex() const
{
    /// I don't know why, but std::max_element returns wrong numbers sometimes!
    //return std::max_element(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(), SymmetrizationIndices_.value_comp())->second;
    unsigned max(0);
    for (auto& kv : SymmetrizationIndices_)
        if (kv.second > max)
            max = kv.second;

    return max;
}

//-------------------------
std::vector<std::shared_ptr<ParticleCombination> > DataAccessor::particleCombinations() const
{
    std::vector<std::shared_ptr<ParticleCombination> > retVal;
    for (auto& kv : SymmetrizationIndices_)
        retVal.push_back(kv.first);

    return retVal;
}

//-------------------------
bool DataAccessor::consistent() const
{
    if (SymmetrizationIndices_.empty()) {
        LOG(ERROR) << "DataAccessor::consistent() - SymmetrizationIndices_ is empty.";
        return false;
    }

    bool result = true;

    for (auto& kv : SymmetrizationIndices_)
        result &= kv.first->consistent();

    return result;
}

//-------------------------
void DataAccessor::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c)
{
    if (SymmetrizationIndices_.find(c) != SymmetrizationIndices_.end())
        // c is already in map
        return;

    // check to see if new member equates to existing member
    for (auto& kv : SymmetrizationIndices_)
        if ((*Equiv_)(kv.first, c)) {
            // equating member found; set index; return
            SymmetrizationIndices_[c] = kv.second;
            return;
        }

    // else assign to current size = highest current index + 1
    SymmetrizationIndices_[c] = maxSymmetrizationIndex() + 1;
}

//-------------------------
void DataAccessor::clearSymmetrizationIndices()
{
    SymmetrizationIndices_.clear();
}

//-------------------------
std::vector<double>& DataAccessor::data(DataPoint& d, unsigned i) const
{
    return d.Data_.at(index()).at(i);
}

//-------------------------
const std::vector<double>& DataAccessor::data(const DataPoint& d, unsigned i) const
{
    return d.Data_.at(index()).at(i);
}

//-------------------------
CalculationStatus& DataAccessor::CalculationStatuses(DataPoint& d, unsigned i)
{
    return d.CalculationStatuses_.at(index()).at(i);
}

//-------------------------
CalculationStatus DataAccessor::CalculationStatuses(DataPoint& d, unsigned i) const
{
    return d.CalculationStatuses_.at(index()).at(i);
}

//-------------------------
InitialStateParticle* DataAccessor::initialStateParticle() const
{
    return InitialStateParticle_;
}


}

