#ifndef __BAT__BAT_YAP_BASE__H
#define __BAT__BAT_YAP_BASE__H

#include <BAT/BCModel.h>

#include <fwd/DecayingParticle.h>
#include <fwd/Model.h>
#include <MassAxes.h>

#include <memory>
#include <numeric>
#include <string>
#include <vector>

class bat_yap_base : public BCModel
{

public:

    bat_yap_base(std::string name, std::unique_ptr<yap::Model> M);

    virtual void MCMCUserInitialize() override
    { LikelihoodCalls_.assign(GetNChains(), 0); }

    unsigned likelihoodCalls() const
    { return std::accumulate(LikelihoodCalls_.begin(), LikelihoodCalls_.end(), 0); }

protected:

    std::shared_ptr<yap::DecayingParticle> isp()
    { return ISP_; }

    std::unique_ptr<yap::Model>& model()
    { return Model_; }

    yap::MassAxes& axes()
    { return Axes_; }

    void increaseLikelihoodCalls(unsigned c)
    { ++LikelihoodCalls_[c]; }

    void increaseLikelihoodCalls()
    { increaseLikelihoodCalls(GetCurrentChain()); }

private:
    yap::MassAxes Axes_;
    std::unique_ptr<yap::Model> Model_;
    std::shared_ptr<yap::DecayingParticle> ISP_;

    std::vector<unsigned> LikelihoodCalls_;

};

#endif
