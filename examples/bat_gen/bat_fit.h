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
    bat_fit(std::string name, std::unique_ptr<yap::Model> M, TTree& t_pars);

    void addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex = "", std::string units = "");
    void addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex = "", std::string units = "");

    /// log likelihood
    double LogLikelihood(const std::vector<double>&) override;

    /// load data from a TTree into a DataSet
    /// \param data DataSet to load into
    /// \param t_mcmc TTree to load from
    /// \param N max number of data points to (attempt to) load
    /// \param lag Lag to apply to iterations when reading from TTree
    size_t loadData(yap::DataSet& data, TTree& t_mcmc, int N = -1, unsigned lag = 1);

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

    integrator_type& integrator()
    { return Integrator_; }

private:

    /// DataSet to fit the model to
    yap::DataSet FitData_;

    /// Partitioning of FitData_
    yap::DataPartitionVector FitPartitions_;

    /// DataSet to calculate model integral with
    yap::DataSet IntegralData_;

    /// Partitioning of IntegralData_
    yap::DataPartitionVector IntegralPartitions_;

    integrator_type Integrator_;

    yap::ModelIntegral Integral_;

    yap::ParameterVector Parameters_;

};

#endif
