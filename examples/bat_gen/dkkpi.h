// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__DKKPI__H
#define __BAT__DKKPI__H

#include <Model.h>
#include <SpinAmplitudeCache.h>

#include <memory>


// This is a dkkpi header file.
// Model source code is located in file yap_test/dkkpi.cxx

// ---------------------------------------------------------
class dkkpi : public yap::Model
{

public:

    // Constructor
    dkkpi(std::unique_ptr<yap::SpinAmplitudeCache> SAC);

};
// ---------------------------------------------------------

#endif
