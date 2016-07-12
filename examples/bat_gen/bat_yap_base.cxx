// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_yap_base.h"

#include <DecayingParticle.h>
#include <Exceptions.h>
#include <Model.h>
#include <SpinAmplitudeCache.h>

// -----------------------
bat_yap_base::bat_yap_base(std::string name, std::unique_ptr<yap::Model> M)
    : BCModel(name),
      Model_(std::move(M))
{
    if (!Model_)
        throw yap::exceptions::Exception("Model nullptr", "bat_yap_base::bat_yap_base");

    Model_->lock();

    if (!Model_->consistent())
        throw yap::exceptions::Exception("Model inconsistent", "bat_yap_base::bat_yap_base");

    auto isps = full_final_state_isp(*Model_);
    if (isps.empty())
        throw yap::exceptions::Exception("no full-final-state initial-state particle in model", "bat_yap_base::bat_yap_base");

    ISP_ = isps[0];
}
