#include "bat_fit.h"

#include "ConstantPrior.h"

#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <SpinAmplitude.h>
#include <VariableStatus.h>

#include <BAT/BCPrior.h>

#include <TTree.h>

void unambiguous_importance_sampler_calculate(yap::ModelIntegral& M, bat_fit::Generator G, unsigned N, unsigned n, unsigned t)
{ yap::ImportanceSampler::calculate(M, G, N, n, t); }

// -----------------------
bat_fit::bat_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : bat_yap_base(name, std::move(M)),
      FitData_(model()->createDataSet()),
      FitPartitions_(1, &FitData_),
      NIntegrationPoints_(0),
      NIntegrationPointsBatchSize_(0),
      NIntegrationThreads_(1),
      Integrator_(integrator_type(unambiguous_importance_sampler_calculate)),
      Integral_(*model()),
      FirstParameter_(-1),
      FirstObservable_(-1)
{
    // create mass axes
    axes() = model()->massAxes(pcs);

    // loop over all free amplitudes
    for (const auto& fa : free_amplitudes(*model())) {
        // ignore fixed free amplitudes
        if (fa->variableStatus() == yap::VariableStatus::fixed)
            continue;

        auto fa_name = to_string(*fa->decayChannel())
            + " L = " + std::to_string(fa->spinAmplitude()->L())
            + " S = " + yap::spin_to_string(fa->spinAmplitude()->twoS())
            + " M = " + yap::spin_to_string(fa->twoM());

        // add real parameter
        AddParameter("real(" + fa_name + ")", -2 * abs(fa->value()), 2 * abs(fa->value()));
        // add imag parameter
        AddParameter("imag(" + fa_name + ")", -2 * abs(fa->value()), 2 * abs(fa->value()));

        AbsPriors_.push_back(new ConstantPrior(0, 2 * abs(fa->value())));
        ArgPriors_.push_back(new ConstantPrior(-180, 180));

        // add amplitude observable
        AddObservable("amp(" + fa_name + ")", 0, 2 * abs(fa->value()));
        // add phase observable
        AddObservable("arg(" + fa_name + ")", -180, 180);

        // add free amplitude to list
        FreeAmplitudes_.push_back(fa);
    }
    
    // // add observables for all fit fractions
    // int N = std::accumulate(Integral_.integrals().begin(), Integral_.integrals().end(), 0,
    //                         [](int n, const yap::IntegralMap::value_type& v)
    //                         {return n + v.second.decayTrees().size();});
    // if (N > 1) {
    unsigned i(0);
    for (const auto& b2_dtvi : Integral_.integrals())
        for (const auto& dt : b2_dtvi.second.decayTrees()) {
            DecayTrees_.push_back(dt);
            AddObservable("fit_frac(" + to_string(*dt->freeAmplitude()) + ")" + std::to_string(i++), 0, 1.1);
        }
    // }

    FirstParameter_ = GetParameters().Size();
    FirstObservable_ = GetObservables().Size();
}

//-------------------------
std::vector<std::vector<unsigned> > find_mass_axes(TTree& t_pars)
{
    std::vector<std::vector<unsigned> > pcs;
    pcs.reserve(t_pars.GetEntries());

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
                                             "find_mass_axes");
        if (parname.size() != 5)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" is not right length "
                                             + "(" + std::to_string(parname.size()) + " != 5)",
                                             "find_mass_axes");
        // read indices:
        std::vector<unsigned> indices;
        indices.reserve(2);
        std::transform(parname.begin() + parname.rfind("_") + 1, parname.end(), std::back_inserter(indices), [](char c) {return std::atoi(&c);});
        // for (size_t i = parname.rfind("_") + 1; i < parname.size(); ++i)
        //     indices.push_back(std::stoi(parname.substr(i, 1)));
        pcs.push_back(indices);
    }

    return pcs;
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
size_t load_data(yap::DataSet& data, const yap::Model& M, const yap::MassAxes& A, double initial_mass, TTree& t_mcmc, int N, unsigned lag)
{
    if (A.empty())
        throw yap::exceptions::Exception("mass axes empty", "load_data");

    // set branch addresses
    std::vector<double> m2;
    m2.reserve(A.size());
    std::for_each(A.begin(), A.end(), std::bind(set_address, std::placeholders::_1, std::ref(m2), std::ref(t_mcmc)));
    if (m2.size() != A.size())
        throw yap::exceptions::Exception("not all mass axes loaded from TTree", "load_data");

    //
    // load data
    //
    int Phase = -1;
    t_mcmc.SetBranchAddress("Phase", &Phase);
    unsigned Iteration;
    t_mcmc.SetBranchAddress("Iteration", &Iteration);

    int n_attempted = 0;
    size_t old_size = data.size();

    for (long long n = 0; n < t_mcmc.GetEntries() and (N < 0 or n_attempted < N); ++n) {
        t_mcmc.GetEntry(n);

        if (Phase <= 0)
            continue;

        if (Iteration % lag != 0)
            continue;

        // if (fabs(m2[0] - 1.35 * 1.35) > 0.1 or m2[1] > 1.55 or m2[1] < 0.58)
        //     continue;

        ++n_attempted;

        auto P = calculate_four_momenta(initial_mass, M, A, m2);
        if (P.empty())
            std::cout << "point is out of phase space!";
        data.push_back(P);
    }

    if (data.empty())
        LOG(INFO) << "No data loaded.";
    else {
        LOG(INFO) << "Loaded " << data.size() - old_size << " data points (" << ((data.size() - old_size) * data[0].bytes() * 1.e-6) << " MB)";
        if (old_size != 0)
            LOG(INFO) << "Total data size now " << data.size() << " points (" << (data.bytes() * 1.e-6) << " MB)";
    }
    return data.size() - old_size;
}

//-------------------------
void bat_fit::setParameters(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i)
        *FreeAmplitudes_[i] = std::complex<double>(p[i * 2], p[i * 2 + 1]);

    yap::set_values(Parameters_.begin(), Parameters_.end(), p.begin() + FirstParameter_, p.end());

    Integrator_(Integral_, IntegrationPointGenerator_, NIntegrationPoints_, NIntegrationPointsBatchSize_, NIntegrationThreads_);

    // calculate fit fractions
    // if (!CalculatedFitFractions_.empty()) {
    unsigned c = GetCurrentChain();
    size_t i = 0;
    for (const auto& b2_dtvi : Integral_.integrals()) {
        auto ff = fit_fractions(b2_dtvi.second);
        for (const auto& f : ff)
            CalculatedFitFractions_[c][i++] = f;
    }
    // }
}

// ---------------------------------------------------------
double bat_fit::LogLikelihood(const std::vector<double>& p)
{
    setParameters(p);
    double L = sum_of_log_intensity(*model(), FitPartitions_, log(integral(Integral_).value()));
    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();
    return L;
}

//-------------------------
double bat_fit::LogAPrioriProbability(const std::vector<double>& p)
{
    double logP = 0;
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i) {
        auto A = std::complex<double>(p[i * 2 + 0], p[i * 2 + 1]);
        logP += AbsPriors_[i]->GetLogPrior(abs(A))
            + ArgPriors_[i]->GetLogPrior(yap::deg(arg(A)))
            - log(abs(A));      // jacobian
    }
    for (size_t i = FreeAmplitudes_.size() * 2; i < GetParameters().Size(); ++i)
        if (GetParameter(i).GetPrior())
            logP += GetParameter(i).GetLogPrior(p[i]);
    return logP;
}

//-------------------------
void bat_fit::CalculateObservables(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i) {
        auto A = std::complex<double>(p[i * 2 + 0], p[i * 2 + 1]);
        GetObservables()[i * 2 + 0] = abs(A);
        GetObservables()[i * 2 + 1] = yap::deg(arg(A));
    }
    unsigned c = GetCurrentChain();
    for (size_t i = 0; i < CalculatedFitFractions_[c].size(); ++i)
        GetObservables()[FreeAmplitudes_.size() * 2 + i] = CalculatedFitFractions_[c][i].value();
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    AddParameter(name + "_re", real(low), real(high), latex.empty() ? latex : "Re(" + latex + ")", units);
    AddParameter(name + "_im", imag(low), imag(high), latex.empty() ? latex : "Im(" + latex + ")", units);
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    AddParameter(name, low, high, latex, units);
}

//-------------------------
size_t bat_fit::findFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A) const
{
    auto it = std::find(FreeAmplitudes_.begin(), FreeAmplitudes_.end(), A);

    if (it == FreeAmplitudes_.end()) {
        LOG(INFO) << "did not find amplitude " << A;
        throw yap::exceptions::Exception("FreeAmplitude not found", "setPrior");
    }

    return (it - FreeAmplitudes_.begin()) * 2;
}

//-------------------------
void bat_fit::setPriors(std::shared_ptr<yap::FreeAmplitude> fa, BCPrior* amp_prior, BCPrior* arg_prior)
{
    if (!amp_prior)
        throw yap::exceptions::Exception("amp_prior is null", "bat_fit::setPrior");
    if (!arg_prior)
        throw yap::exceptions::Exception("phase_prior is null", "bat_fit::setPrior");
    
    auto i = findFreeAmplitude(fa);

    AbsPriors_[i / 2] = amp_prior;
    ArgPriors_[i / 2] = arg_prior;
}

//-------------------------
void bat_fit::setRealImagRanges(std::shared_ptr<yap::FreeAmplitude> fa, double real_low, double real_high, double imag_low, double imag_high)
{
    auto i = findFreeAmplitude(fa);
    GetParameter(i).SetLimits(real_low, real_high);
    GetParameter(i + 1).SetLimits(imag_low, imag_high);
}
    
//-------------------------
void bat_fit::setAbsArgRanges(std::shared_ptr<yap::FreeAmplitude> fa, double abs_low, double abs_high, double arg_low, double arg_high)
{
    auto i = findFreeAmplitude(fa);
    GetObservable(i).SetLimits(abs_low, abs_high);
    GetObservable(i + 1).SetLimits(arg_low, arg_high);
}

//-------------------------
void bat_fit::fix(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase)
{
    auto i = findFreeAmplitude(A);
    auto a = std::polar(amp, yap::rad(phase));
    GetParameter(i).SetLimits(real(a) - 1, real(a) + 1);
    GetParameter(i).Fix(real(a));
    GetParameter(i + 1).SetLimits(imag(a) - 1, imag(a) + 1);
    GetParameter(i + 1).Fix(imag(a));
    GetObservable(i).SetLimits(std::max<double>(0, amp - 1), amp + 1);
    GetObservable(i + 1).SetLimits(phase - 15, phase + 15);
}

//-------------------------
void bat_fit::MCMCUserInitialize()
{
    bat_yap_base::MCMCUserInitialize();
    // if (!CalculatedFitFractions_.empty())
    CalculatedFitFractions_.assign(GetNChains(), yap::RealIntegralElementVector(DecayTrees_.size()));
}

