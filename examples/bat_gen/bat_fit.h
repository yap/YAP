#ifndef __BAT__BAT_FIT__H
#define __BAT__BAT_FIT__H

#include "bat_yap_base.h"

#include <fwd/Model.h>

#include <DataSet.h>

#include <memory>
#include <string>
#include <vector>

class TTree;

class bat_fit : public bat_yap_base
{

public:

    bat_fit(std::string name, std::unique_ptr<yap::Model> M,
            TTree& t_mcmc, TTree& t_pars, int N = -1, unsigned lag = 1);

    double LogLikelihood(const std::vector<double>&) override;

private:

    yap::DataSet Data_;

};

#endif
