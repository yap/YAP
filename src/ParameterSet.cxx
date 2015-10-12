#include "ParameterSet.h"

#include "logging.h"

namespace yap {

//-------------------------
ParameterSet& ParameterSet::operator=(std::initializer_list<double> pars)
{
    //Parameters_ = pars; synchronizeVariableStatuses(); return *this;

    std::vector<double> parsVec(pars);

    if (parsVec.size() != Parameters_.size())
        LOG(FATAL) << "ParameterSet& operator= : pars have wrong size";

    for (unsigned i = 0; i < parsVec.size(); ++i) {
        if (parsVec[i] != Parameters_[i]) {
            Parameters_[i] = parsVec[i];
            VariableStatuses_[i] = kChanged;
        }
    }

    return *this;
}

}
