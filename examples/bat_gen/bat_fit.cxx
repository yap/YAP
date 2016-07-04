#include "bat_fit.h"

#include <DataSet.h>
#include <DecayingParticle.h>
#include <FourVector.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>

#include <TTree.h>

// -----------------------
bat_fit::bat_fit(std::string name, std::unique_ptr<yap::Model> M, TTree& t_mcmc, TTree& t_pars,
                 int N, unsigned lag)
    : bat_yap_base(name, std::move(M)),
      Data_(model()->createDataSet())
{
    unsigned n_fsp = model()->finalStateParticles().size();
    unsigned n_dof = 3 * n_fsp - 7;

    //
    // find mass axes
    // and set t_mcmc parameter branch addresses
    //
    std::vector<std::vector<unsigned> > pcs;
    pcs.reserve(n_dof);
    std::vector<double> m2;
    m2.reserve(n_dof);

    bool parameter;
    char c_parname[10];
    t_pars.SetBranchAddress("parameter", &parameter);
    t_pars.SetBranchAddress("name", &c_parname);
    for (unsigned n = 0; n < t_pars.GetEntries(); ++n) {
        t_pars.GetEntry(n);
        if (!parameter)
            continue;
        std::string parname(c_parname);
        // parameter names should be m2_ij
        if (parname.find("m2_") != 0)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" does not match \"m2_ij\"",
                                             "bat_fit::bat_fit");
        if (parname.size() != 5)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" is not right length "
                                             + "(" + std::to_string(parname.size()) + " != 5)",
                                             "bat_fit::bat_fit");
        // read indices:
        std::vector<unsigned> indices;
        indices.reserve(2);
        for (size_t i = parname.rfind("_") + 1; i < parname.size(); ++i) {
            indices.push_back(std::stoi(parname.substr(i, 1)));
            if (indices.back() >= n_fsp)
                throw yap::exceptions::Exception("index out of range (" + std::to_string(indices.back()) + " >= "
                                                 + std::to_string(n_fsp) + ")",
                                                 "bat_fit::bat_fit");
        }
        pcs.push_back(indices);

        // set branch address
        m2.push_back(0);
        t_mcmc.SetBranchAddress(parname.data(), &m2.back());
    }

    if (pcs.size() != n_dof or m2.size() != n_dof)
        throw yap::exceptions::Exception("insufficient axes found "
                                         "(" + std::to_string(pcs.size()) + " < " + std::to_string(n_dof) + ")",
                                         "bat_fit::bat_fit");

    // create mass axes
    axes() = model()->massAxes(pcs);

    //
    // load data
    //
    int Phase = -1;
    t_mcmc.SetBranchAddress("Phase", &Phase);
    unsigned Iteration;
    t_mcmc.SetBranchAddress("Iteration", &Iteration);

    long long n = 0;
    // find first entry of main run
    while (Phase <= 0 and n < t_mcmc.GetEntries())
        t_mcmc.GetEntry(n++);

    auto nN = (N < 0) ? t_mcmc.GetEntries() : std::min(n + N, t_mcmc.GetEntries());

    while (n < nN) {
        t_mcmc.GetEntry(n++);

        if (Phase <= 0)
            continue;

        if (Iteration % lag != 0)
            continue;

        auto P = model()->calculateFourMomenta(axes(), m2, isp()->mass()->value());
        if (P.empty())
            std::cout << "point is out of phase space!";
        Data_.push_back(P);
    }

    LOG(INFO) << "Data loaded with " << Data_.size() << " points";
    if (!Data_.points().empty())
        LOG(INFO) << "Total data size = " << (Data_.size() * Data_[0].bytes() * 1.e-6) << " MB";
}

// ---------------------------------------------------------
double bat_fit::LogLikelihood(const std::vector<double>&)
{
    unsigned c = GetCurrentChain();
    double L = sum_of_log_intensity(*model(), Data_);
    // model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls(c);
    return L;
}
