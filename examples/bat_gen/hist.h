#ifndef __hist__h
#define __hist__h

#include "DataPoint.h"
#include "DataSet.h"
#include "FourMomenta.h"
#include "MassAxes.h"
#include "MassRange.h"
#include "ParticleCombination.h"

#include <TH1D.h>
#include <TH2D.h>

#include <iostream>

TH2D* hist2(const yap::FourMomenta& P, const yap::MassAxes& A, const std::vector<yap::MassRange> m2r, const yap::DataSet& data)
{
    if (A.size() != 2 or m2r.size() != 2)
        throw;

    int nbins = sqrt(data.size() / 10.);

    TH2D* h2 = new TH2D("h2", (";" + yap::indices_string(*A[0]) + ";" + yap::indices_string(*A[1])).data(),
                        nbins, m2r[0][0], m2r[0][1], nbins, m2r[1][0], m2r[1][1]);
    for (const auto& d : data)
        h2->Fill(P.m2(d, A[0]), P.m2(d, A[1]));
    return h2;
}

TH1D* hist1(const yap::FourMomenta& P, const yap::MassAxes::value_type& a, const yap::MassRange m2r, const yap::DataSet& data)
{
    int nbins = data.size() / 100.;

    TH1D* h1 = new TH1D("h1", (";" + yap::indices_string(*a)).data(),
                        nbins, m2r[0], m2r[1]);
    for (const auto& d : data) {
        h1->Fill(P.m2(d, a));
    }
    return h1;
}

#endif
