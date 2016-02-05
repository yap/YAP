/*
 *  Created by Phil on 28/04/2011.
 *  Copyright 2010 Two Blue Cubes Ltd. All rights reserved.
 *
 *  Distributed under the Boost Software License, Version 1.0. (See accompanying
 *  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef TWOBLUECUBES_CATCH_APPROX_HPP_INCLUDED
#define TWOBLUECUBES_CATCH_APPROX_HPP_INCLUDED

#include "catch_tostring.h"

#include <cmath>
#include <limits>

namespace Catch {
namespace Detail {

    template <typename T>
    class Approx {
    public:
        explicit Approx ( T value)
        :   m_epsilon( std::numeric_limits<float>::epsilon()*100 ),
            m_scale( 1.0 ),
            m_value( value )
        {}

        Approx( Approx const& other )
        :   m_epsilon( other.m_epsilon ),
            m_scale( other.m_scale ),
            m_value( other.m_value )
        {}

        static Approx custom() {
            return Approx( 0 );
        }

        Approx operator()( T value ) {
            Approx approx( value );
            approx.epsilon( m_epsilon );
            approx.scale( m_scale );
            return approx;
        }

        friend bool operator == ( const T& lhs, Approx const& rhs ) {
            return operator==( rhs, lhs ); 
        }

        friend bool operator != ( const T& lhs, Approx const& rhs ) {
            return !operator==( lhs, rhs );
        }

        friend bool operator != ( Approx const& lhs, const T& rhs ) {
            return !operator==( rhs, lhs );
        }

        Approx& epsilon( double newEpsilon ) {
            m_epsilon = newEpsilon;
            return *this;
        }

        Approx& scale( double newScale ) {
            m_scale = newScale;
            return *this;
        }

        std::string toString() const {
            std::ostringstream oss;
            oss << "Approx( " << Catch::toString( m_value ) << " )";
            return oss.str();
        }

        bool check (const double& lhs, const double& rhs) const {
                return fabs( lhs - rhs ) < m_epsilon * (m_scale + (std::max)( fabs(lhs), fabs(rhs) ) );
        }

        friend operator == (const T&, const Approx<T>&);

    private:
        double m_epsilon;
        double m_scale;
        T m_value;
    };
}

bool operator == ( double lhs, Approx<double> const& rhs ) {
    return rhs.check(lhs, rhs.m_value);
}

bool operator == ( std::complex<double> lhs, Approx<std::complex<double> > const& rhs ) {
    return rhs.check(lhs[0], rhs.m_value[0]) and rhs.check(lhs[1], rhs.m_value[1]);
}


template<>
inline std::string toString<Detail::Approx>( Detail::Approx const& value ) {
    return value.toString();
}

} // end namespace Catch

#endif // TWOBLUECUBES_CATCH_APPROX_HPP_INCLUDED
