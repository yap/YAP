// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI__H
#define __BAT__D3PI__H

#include <Model.h>
#include <SpinAmplitudeCache.h>

#include <memory>

std::unique_ptr<yap::Model> d3pi(std::unique_ptr<yap::SpinAmplitudeCache> SAC);

#endif
