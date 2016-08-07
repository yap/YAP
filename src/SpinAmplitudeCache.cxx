#include "SpinAmplitudeCache.h"

#include "Exceptions.h"

namespace yap {

//-------------------------
SpinAmplitudeCache::weak_ptr_type SpinAmplitudeCache::find(unsigned two_J, const SpinVector& two_j, unsigned L, unsigned two_S) const
{
    for (const auto& wp : *this) {
        if (wp.expired())
            continue;
        auto sp = wp.lock();
        if (sp->initialTwoJ() == two_J and sp->finalTwoJ() == two_j and sp->L() == L and sp->twoS() == two_S)
            return wp;
    }
    return weak_ptr_type();
}

//-------------------------
std::shared_ptr<SpinAmplitude> SpinAmplitudeCache::spinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned L, unsigned two_S)
{
    // search for existing spin amplitude
    auto wp = find(two_J, two_j, L, two_S);
    // if valid one found, return it
    if (!wp.expired())
        return wp.lock();
    // else make a new one
    if (!Model_)
        throw exceptions::Exception("Model is nullptr", "SpinAmplitudeCache::spinAmplitude");
    // if more than two-body
    if (two_j.size() > 2)
        return operator[](unit(*Model_, two_J, two_j, L, two_S));
    // else
    return operator[](create(*Model_, two_J, two_j, L, two_S));
}

//-------------------------
bool SpinAmplitudeCache::consistent() const
{
    bool C = true;
    for (auto it = begin(); it != end(); ++it)
        if (!it->expired())
            C &= it->lock()->consistent();
    return C;
}

//-------------------------
void SpinAmplitudeCache::setModel(Model& model)
{
    if (Model_ != nullptr)
        throw exceptions::Exception("Model already set.", "SpinAmplitudeCache::setModel");
    Model_ = &model;
}

}
