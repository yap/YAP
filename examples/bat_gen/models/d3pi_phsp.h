// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI_PHSP__H
#define __BAT__D3PI_PHSP__H

#include "find_pdl_file.h"

#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Model.h>
#include <ParticleTable.h>
#include <PDL.h>

#include <memory>

inline std::unique_ptr<yap::Model> d3pi_phsp(std::unique_ptr<yap::Model> M)
{
    auto T = yap::read_pdl_file(find_pdl_file());

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    M->setFinalState(piPlus, piMinus, piPlus);

    auto D = DecayingParticle::create(T["D+"], 3 /*Gev^-1*/);
    D->addWeakDecay(piPlus, piMinus, piPlus);

    return M;
}

#endif
