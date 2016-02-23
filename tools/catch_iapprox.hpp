#ifndef yap_catch_IApprox_h
#define yap_catch_IApprox_h

#include "catch_vapprox.hpp"

#include <stdexcept>

namespace Catch {
namespace Detail {

    template <class T>
    class IApprox : public VApprox<T>
    {
    public:
        IApprox(const T& v, double eps = std::numeric_limits<float>::epsilon() * 100) :
            VApprox<T>(v, eps)
            {}

        bool operator==(const T& rhs) const override
            {
                bool R = true;
                auto it = this->value().begin();
                auto rhs_it = rhs.begin();
                while (it != this->value().end() and rhs_it != rhs.end()) {
                    R &= this->equal((double)*it, (double)*rhs_it);
                    ++it;
                    ++rhs_it;
                }
                if (it != this->value().end() or rhs_it != rhs.end())
                    throw std::exception();
                return R;
            }

    };

}
}

#endif
