#include "logging.h"
INITIALIZE_EASYLOGGINGPP

// #include "DataSet.h"
// #include "DataPoint.h"

#include "ParticleCombination.h"
#include "ParticleIndex.h"

#include <iostream>
#include <memory>
#include <vector>

// #include <TLorentzVector.h>

int main( int argc, char** argv)
{

    // FSP's
    std::shared_ptr<yap::ParticleCombination> pi1 = std::make_shared<yap::ParticleCombination>(0);
    std::shared_ptr<yap::ParticleCombination> pi2 = std::make_shared<yap::ParticleCombination>(1);
    std::shared_ptr<yap::ParticleCombination> pi3 = std::make_shared<yap::ParticleCombination>(2);
    std::shared_ptr<yap::ParticleCombination> pi4 = std::make_shared<yap::ParticleCombination>(3);

    std::shared_ptr<yap::ParticleCombination> rho1 = std::make_shared<yap::ParticleCombination>();
    rho1->addDaughter(pi1);
    rho1->addDaughter(pi2);

    std::shared_ptr<yap::ParticleCombination> rho2 = std::make_shared<yap::ParticleCombination>();
    rho2->addDaughter(pi3);
    rho2->addDaughter(pi4);

    std::shared_ptr<yap::ParticleCombination> D1 = std::make_shared<yap::ParticleCombination>();
    D1->addDaughter(rho1);
    D1->addDaughter(rho2);

    std::vector<yap::ParticleIndex> I = D1->indices();

    for (unsigned i = 0; i < I.size(); ++i)
        std::cout << i << "\t" << static_cast<int>(I[i]) << std::endl;

    // std::vector<TLorentzVector> P;
    // P.push_back(TLorentzVector(0, 1, 2, 3));
    // P.push_back(TLorentzVector(1, 2, 3, 0));
    // P.push_back(TLorentzVector(2, 3, 0, 1));
    // P.push_back(TLorentzVector(3, 0, 1, 2));

    // yap::DataSet DS;
    // DS.addDataPoint(yap::DataPoint(P)); // move

    // yap::DataPoint D(P);
    // DS.addDataPoint(D);         // copy and move
    // DS.addDataPoint(std::move(D)); // move

    // P.pop_back();
    // DS.addDataPoint(yap::DataPoint(P)); // fail



}
