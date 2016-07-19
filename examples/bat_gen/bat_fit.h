#ifndef __BAT__BAT_FIT__H
#define __BAT__BAT_FIT__H

#include "bat_yap_base.h"

#include <fwd/DataPartition.h>
#include <fwd/Model.h>
#include <fwd/Parameter.h>

#include <DataSet.h>
#include <ModelIntegral.h>

#include <memory>
#include <string>
#include <vector>

class TTree;

class bat_fit : public bat_yap_base
{

public:

    /// constructor
    bat_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs);

    /// add a complex parameter from the model to bat
    void addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex = "", std::string units = "");

    /// ada a real parameter from the model to bat
    void addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex = "", std::string units = "");

    /// log likelihood
    double LogLikelihood(const std::vector<double>&) override;

    /// \return FitData_
    yap::DataSet& fitData()
    { return FitData_; }

    /// \return FitPartitions_
    yap::DataPartitionVector& fitPartitions()
    { return FitPartitions_; }

    /// \return NormalizationData_
    yap::DataSet& integralData()
    { return IntegralData_; }

    /// \return NormalizationData_
    yap::DataPartitionVector& integralPartitions()
    { return IntegralPartitions_; }

    /// \typedef integrator_type
    /// convienence typedef
    using integrator_type = std::function<void(yap::ModelIntegral&, yap::DataPartitionVector&)>;

    /// \return the integrator
    integrator_type& integrator()
    { return Integrator_; }

protected:

    /// DataSet to fit the model to
    yap::DataSet FitData_;

    /// Partitioning of FitData_
    yap::DataPartitionVector FitPartitions_;

    /// DataSet to calculate model integral with
    yap::DataSet IntegralData_;

    /// Partitioning of IntegralData_
    yap::DataPartitionVector IntegralPartitions_;

    /// function to integrate using
    integrator_type Integrator_;

    /// stores integral result
    yap::ModelIntegral Integral_;

    /// vector of parameters to set in model
    yap::ParameterVector Parameters_;

};

/// load data from a TTree into a DataSet
/// \param data DataSet to load into
/// \param M Model to load with
/// \param A MassAxes to load with
/// \param t_mcmc TTree to load from
/// \param N max number of data points to (attempt to) load
/// \param lag Lag to apply to iterations when reading from TTree
size_t load_data(yap::DataSet& data, const yap::Model& M,
                 const yap::MassAxes& A, double initial_mass, TTree& t_mcmc,
                 int N = -1, unsigned lag = 1);

/// find mass axes from TTree of parameters
std::vector<std::vector<unsigned> > find_mass_axes(TTree& t_pars);

#endif
