#ifndef __BAT__BAT_FIT__H
#define __BAT__BAT_FIT__H

#include "bat_yap_base.h"

#include <fwd/DataPartition.h>
#include <fwd/FourVector.h>
#include <fwd/FreeAmplitude.h>
#include <fwd/IntegralElement.h>
#include <fwd/Model.h>
#include <fwd/Parameter.h>

#include <DataSet.h>
#include <ModelIntegral.h>

#include <memory>
#include <functional>
#include <string>
#include <random>
#include <vector>

class TTree;

class bat_fit : public bat_yap_base
{

public:

    /// constructor
    /// \param name name of bat model
    /// \param M yap::model
    /// \param pcs vector<vector<unsigned>> defining mass axes
    bat_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs = {});

    /// add a complex parameter from the model to bat
    void addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex = "", std::string units = "");

    /// ada a real parameter from the model to bat
    void addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex = "", std::string units = "");

    /// set the range and prior for a FreeAmplitude
    void setPrior(std::shared_ptr<yap::FreeAmplitude> A, BCPrior* amp_prior, BCPrior* phase_prior);

    /// set the range and a constant prior for a FreeAmplitude
    void setPrior(std::shared_ptr<yap::FreeAmplitude> A, double amp_low, double amp_high, double phase_low, double phase_high);

    /// fix a FreeAmplitude
    void fix(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase);

    /// log likelihood
    double LogLikelihood(const std::vector<double>& p) override;

    /// calculate  observables
    void CalculateObservables(const std::vector<double>& p) override;

    /// \return FitData_
    yap::DataSet& fitData()
    { return FitData_; }

    /// \return FitPartitions_
    yap::DataPartitionVector& fitPartitions()
    { return FitPartitions_; }

    /// set parameters of integration
    /// \param N number of points to use for integration
    /// \param n batch size for integration
    void setNIntegrationPoints(unsigned N, unsigned n)
    { NIntegrationPoints_ = N; NIntegrationPointsBatchSize_ = n; }

    /// \typedef Generator
    /// function for generating new points for integration
    using Generator = std::function<std::vector<yap::FourVector<double> >()>;

    /// \return IntegrationPointGenerator_
    Generator& integrationPointGenerator()
    { return IntegrationPointGenerator_; }

    /// \typedef integrator_type
    /// convienence typedef
    using integrator_type = std::function<void(yap::ModelIntegral&, Generator, unsigned, unsigned)>;

    /// \return the integrator
    integrator_type& integrator()
    { return Integrator_; }

    /// init size of CalculatedFitFractions_
    void MCMCUserInitialize() override;

    /// set parameters into model
    void setParameters(const std::vector<double>& p);

    /// find the position in the parameter list of the first element of a free amplitude
    size_t findFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A) const;

protected:

    /// DataSet to fit the model to
    yap::DataSet FitData_;

    /// Partitioning of FitData_
    yap::DataPartitionVector FitPartitions_;

    /// Number of points to integrate with
    unsigned NIntegrationPoints_;

    /// Batch size for generating integration points
    unsigned NIntegrationPointsBatchSize_;

    /// generator for integration
    Generator IntegrationPointGenerator_;

    /// function to integrate using
    integrator_type Integrator_;

    /// stores integral result
    yap::ModelIntegral Integral_;

    /// vector of parameters to set in model
    yap::ParameterVector Parameters_;

    /// offset of where first user-set parameter is
    int FirstParameter_;

    /// list of decay trees integrated over
    yap::DecayTreeVector DecayTrees_;

    /// Free amplitudes of model to set
    yap::FreeAmplitudeVector FreeAmplitudes_;

    /// Calculated fit fractions (for observables)
    std::vector<yap::RealIntegralElementVector> CalculatedFitFractions_;

};

/// load data from a TTree into a DataSet
/// \param data DataSet to load into
/// \param M Model to load with
/// \param A MassAxes to load with
/// \param t_mcmc TTree to load from
/// \param N max number of data points to (attempt to) load
/// \param lag Lag to apply to iterations when reading from TTree
/// \param eps Amount to smear momenta by
size_t load_data(yap::DataSet& data, const yap::Model& M,
                 const yap::MassAxes& A, double initial_mass, TTree& t_mcmc,
                 int N = -1, unsigned lag = 1);

/// find mass axes from TTree of parameters
std::vector<std::vector<unsigned> > find_mass_axes(TTree& t_pars);

#endif
