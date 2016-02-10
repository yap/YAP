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

// This is a d3pi header file.
// Model source code is located in file yap_test/d3pi.cxx

// ---------------------------------------------------------
class d3pi : public yap::Model
{

public:

    // Constructor
    d3pi(std::unique_ptr<yap::SpinAmplitudeCache> SAC);

};
// ---------------------------------------------------------

#endif
