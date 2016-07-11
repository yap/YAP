#include "bat_fit.h"

#include <DataSet.h>
#include <DecayingParticle.h>
#include <FourVector.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>

#include <TTree.h>

void unambiguous_importance_sampler_calculate(yap::ModelIntegral& M, yap::DataPartition& D)
{ yap::ImportanceSampler::calculate(M, D); }

// -----------------------
bat_fit::bat_fit(std::string name, std::unique_ptr<yap::Model> M, TTree& t_pars)
    : bat_yap_base(name, std::move(M)),
      FitData_(model()->createDataSet()),
      NormalizationData_(model()->createDataSet()),
      Integrator_(integrator_type(unambiguous_importance_sampler_calculate))
{
    unsigned n_fsp = model()->finalStateParticles().size();
    unsigned n_dof = 3 * n_fsp - 7;

    //
    // find mass axes
    //
    std::vector<std::vector<unsigned> > pcs;
    pcs.reserve(n_dof);

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

    }

    if (pcs.size() != n_dof)
        throw yap::exceptions::Exception("insufficient axes found "
                                         "(" + std::to_string(pcs.size()) + " < " + std::to_string(n_dof) + ")",
                                         "bat_fit::bat_fit");

    // create mass axes
    axes() = model()->massAxes(pcs);
}

//-------------------------
void set_address(const yap::MassAxes::value_type& a,
                 std::vector<double>& m2,
                 TTree& t_mcmc)
{
    m2.push_back(0);
    t_mcmc.SetBranchAddress(indices_string(*a, "m2_", "").data(), &m2.back());
}

//-------------------------
size_t bat_fit::loadData(yap::DataSet& data, TTree& t_mcmc, int N, unsigned lag)
{
    if (axes().empty())
        throw yap::exceptions::Exception("mass axes empty", "bat_fit::loadData");

    // set branch addresses
    std::vector<double> m2;
    m2.reserve(axes().size());
    std::for_each(axes().begin(), axes().end(), std::bind(set_address, std::placeholders::_1, std::ref(m2), std::ref(t_mcmc)));
    if (m2.size() != axes().size())
        throw yap::exceptions::Exception("not all mass axes loaded from TTree", "bat_fit::loadData");

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

    int n_attempted = 0;
    size_t old_size = data.size();

    for (; n < t_mcmc.GetEntries() and (N < 0 or n_attempted < N); ++n) {
        t_mcmc.GetEntry(n);

        if (Phase <= 0)
            continue;

        if (Iteration % lag != 0)
            continue;

        ++n_attempted;

        auto P = model()->calculateFourMomenta(axes(), m2, isp()->mass()->value());
        if (P.empty())
            std::cout << "point is out of phase space!";
        data.push_back(P);
    }

    if (data.empty())
        LOG(INFO) << "No data loaded.";
    else {
        LOG(INFO) << "Loaded " << data.size() - old_size << " data points (" << ((data.size() - old_size) * data[0].bytes() * 1.e-6) << " MB)";
        if (old_size != 0)
            LOG(INFO) << "Total data size now " << data.size() << " points (" << (data.size() * data[0].bytes() * 1.e-6) << " MB)";
    }
    return data.size() - old_size;
}

// ---------------------------------------------------------
double bat_fit::LogLikelihood(const std::vector<double>&)
{
    double L = sum_of_log_intensity(*model(), FitData_);
    
    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();
    return L;
}
