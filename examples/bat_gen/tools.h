#ifndef __tools_h__
#define __tools_h__

#include <FourVector.h>
#include <Model.h>

#include <TLorentzVector.h>

#include <memory>

template <typename T>
std::unique_ptr<yap::Model> yap_model()
{ return std::make_unique<yap::Model>(std::make_unique<T>()); }

inline const double quad(std::vector<double> S)
{ return sqrt(std::accumulate(S.begin(), S.end(), 0., [](double a, double s) {return a + s * s;})); }

template <typename ... Types>
constexpr double quad(double s0, Types ... additional)
{ return quad({s0, additional...}); }

inline yap::FourVector<double> convert(const TLorentzVector& p)
{ return yap::FourVector<double>({p.E(), p.X(), p.Y(), p.Z()}); }

#endif
