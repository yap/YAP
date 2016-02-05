/*
 *  Created by Phil on 28/04/2011.
 *  Copyright 2010 Two Blue Cubes Ltd. All rights reserved.
 *
 *  Distributed under the Boost Software License, Version 1.0. (See accompanying
 *  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef TWOBLUECUBES_CATCH_CApprox_HPP_INCLUDED
#define TWOBLUECUBES_CATCH_CApprox_HPP_INCLUDED

//#include "catch_tostring.h"

#include <cmath>
#include <complex>
#include <limits>

namespace Catch {
namespace Detail {

    class CApprox {
    public:
        explicit CApprox ( std::complex<double> value )
        :   m_epsilon( std::numeric_limits<float>::epsilon()*100 ),
            m_scale( 1.0 ),
            m_value( value )
        {}
        
        explicit CApprox ( std::complex<double> value, double epsilon )
        :   m_epsilon( epsilon ),
            m_scale( 1.0 ),
            m_value( value )
        {}

        CApprox( CApprox const& other )
        :   m_epsilon( other.m_epsilon ),
            m_scale( other.m_scale ),
            m_value( other.m_value )
        {}

        static CApprox custom() {
            return CApprox( 0 );
        }
        
        CApprox operator()( std::complex<double> value ) {
            CApprox CApprox( value );
            CApprox.epsilon( m_epsilon );
            CApprox.scale( m_scale );
            return CApprox;
        }

            
        friend bool operator == ( std::complex<double> lhs, CApprox const& rhs ) {
            // Thanks to Richard Harris for his help refining this formula
            return fabs( lhs.real() - rhs.m_value.real() ) < rhs.m_epsilon * (rhs.m_scale + (std::max)( fabs(lhs.real()), fabs(rhs.m_value.real()) ) )
                   and fabs( lhs.imag() - rhs.m_value.imag() ) < rhs.m_epsilon * (rhs.m_scale + (std::max)( fabs(lhs.imag()), fabs(rhs.m_value.imag()) ) );
        }

        friend bool operator == ( CApprox const& lhs, std::complex<double> rhs ) {
            return operator==( rhs, lhs );
        }

        friend bool operator != ( std::complex<double> lhs, CApprox const& rhs ) {
            return !operator==( lhs, rhs );
        }

        friend bool operator != ( CApprox const& lhs, std::complex<double> rhs ) {
            return !operator==( rhs, lhs );
        }
        

        CApprox& epsilon( double newEpsilon ) {
            m_epsilon = newEpsilon;
            return *this;
        }

        CApprox& scale( double newScale ) {
            m_scale = newScale;
            return *this;
        }

        std::string toString() const {
            std::ostringstream oss;
            oss << "CApprox( " << m_value << " )";
            return oss.str();
        }

    private:
        double m_epsilon;
        double m_scale;
        std::complex<double> m_value;
    };
}

template<>
inline std::string toString<Detail::CApprox>( Detail::CApprox const& value ) {
    return value.toString();
}

} // end namespace Catch

#endif // TWOBLUECUBES_CATCH_CApprox_HPP_INCLUDED
