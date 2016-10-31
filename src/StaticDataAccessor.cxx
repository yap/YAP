#include "StaticDataAccessor.h"

#include "Model.h"

namespace yap {

//-------------------------
void StaticDataAccessor::registerWithModel()
{
    DataAccessor::registerWithModel();
    // this is in the model's data accessors
    // and not yet in its static data accessors, add it
    if (model()->dataAccessors().find(this) != model()->dataAccessors().end() and
        std::find(model()->StaticDataAccessors_.begin(), model()->StaticDataAccessors_.end(), this) == model()->StaticDataAccessors_.end())
        addToStaticDataAccessors();
}

//-------------------------
void StaticDataAccessor::addToStaticDataAccessors()
{
    const_cast<Model*>(model())->StaticDataAccessors_.push_back(this);
}

//-------------------------
void remove_expired(StaticDataAccessorVector& S)
{
    for (auto it = S.begin(); it != S.end(); )
        if (!*it) it = S.erase(it);
        else ++it;
}

}
