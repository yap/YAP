/*
 *  Extending catch::Approx class which has
 *  Copyright 2010 Two Blue Cubes Ltd. All rights reserved.
 */
#ifndef yap_catch_VApprox_h
#define yap_catch_VApprox_h

//#include "catch_tostring.h"

#include <cmath>
#include <limits>

namespace Catch {
namespace Detail {

    template <typename T>
    class VApprox
    {
    public:

        explicit VApprox(const T& v, double eps = std::numeric_limits<float>::epsilon() * 100) :
            m_value(v), m_epsilon(eps), m_scale( 1.0 )
            {}

        bool equal(const double& lhs, const double& rhs) const
            { return fabs(lhs - rhs) < m_epsilon * (m_scale + (std::max)( fabs(lhs), fabs(rhs))); }

        virtual bool operator==(const T& rhs) const = 0;

        friend bool operator==(const T& lhs, const VApprox& rhs)
            { return rhs.operator==(lhs); }

        friend bool operator!=(const T& lhs, const VApprox& rhs)
            { return !rhs.operator==(lhs); }
        
        friend bool operator!=(const VApprox& lhs, const T& rhs)
            { return !rhs.operator==(rhs); }
        
        VApprox& epsilon(double newEpsilon)
            { m_epsilon = newEpsilon; return *this; }

        VApprox& scale(double newScale)
            { m_scale = newScale; return *this; }

        virtual std::string toString() const
            {
                using std::to_string;
                std::ostringstream oss;
                oss << "CApprox( " << to_string(m_value) << " )";
                return oss.str();
            }

        const T& value() const
            { return m_value; }

    private:
        T m_value;
        double m_epsilon;
        double m_scale;
    };
}

    template <typename T>
    inline std::string toString(const Detail::VApprox<T>& value) {
        return value.toString();
}

} // end namespace Catch

#endif
