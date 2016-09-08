#ifndef __tools_h__
#define __tools_h__

#include <make_unique.h>
#include <Model.h>

#include <memory>

template <typename T>
std::unique_ptr<yap::Model> yap_model(bool include_phsp_factors = false)
{ return std::make_unique<yap::Model>(std::make_unique<T>(), include_phsp_factors); }

const double quad(std::vector<double> S)
{ return sqrt(std::accumulate(S.begin(), S.end(), 0., [](double a, double s) {return a + s * s;})); }

template <typename ... Types>
constexpr double quad(double s0, Types ... additional)
{ return quad({s0, additional...}); }


#endif
