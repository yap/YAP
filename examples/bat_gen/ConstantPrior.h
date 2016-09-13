#ifndef __ConstantPrior__h
#define __ConstantPrior__h

#include <BAT/BCConstantPrior.h>

#include <limits>

class ConstantPrior : public BCConstantPrior
{
public:
    /// constructor
    ConstantPrior(double low, double high)
        : BCConstantPrior(low, high), Low_(low), High_(high) {}

    /// log prior
    virtual double GetLogPrior(double p)
    {
        return (p > Low_ and p < High_) ?
            BCConstantPrior::GetLogPrior(p)
            : 
            -std::numeric_limits<double>::infinity();
    }

private:

    double Low_;
    double High_;

};

#endif
