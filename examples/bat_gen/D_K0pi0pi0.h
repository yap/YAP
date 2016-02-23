// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__DKSPI0PI0__H
#define __BAT__DKSPI0PI0__H

#include <DecayingParticle.h>
#include <Model.h>
#include <SpinAmplitudeCache.h>

#include <memory>

std::unique_ptr<yap::Model> D_K0pi0pi0(std::unique_ptr<yap::SpinAmplitudeCache> SAC);

#endif
