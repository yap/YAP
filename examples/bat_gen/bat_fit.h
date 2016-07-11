#ifndef __BAT__BAT_FIT__H
#define __BAT__BAT_FIT__H

#include "bat_yap_base.h"

#include <fwd/DataPartition.h>
#include <fwd/Model.h>
#include <fwd/ModelIntegral.h>

#include <DataSet.h>

#include <memory>
#include <string>
#include <vector>

class TTree;

class bat_fit : public bat_yap_base
{

public:
    
    /// constructor
    bat_fit(std::string name, std::unique_ptr<yap::Model> M,
            TTree& t_pars);

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

    /// \return NormalizationData_
    yap::DataSet& normalizationData()
    { return NormalizationData_; }

    /// \typedef integrator_type
    /// convienence typedef
    using integrator_type = std::function<void(yap::ModelIntegral&, yap::DataPartition&)>;
    
    integrator_type& integrator()
    { return Integrator_; }

private:

    /// DataSet to fit the model to
    yap::DataSet FitData_;

    /// DataSet to calculate model integral with
    yap::DataSet NormalizationData_;

    integrator_type Integrator_;

};

#endif
