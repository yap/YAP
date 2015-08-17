#include "logging.h"
INITIALIZE_EASYLOGGINGPP

#include "ParticleCombination.h"
#include "ParticleIndex.h"

#include <iostream>
#include <memory>
#include <set>
#include <vector>


int main( int argc, char** argv)
{

    // FSP's
    std::shared_ptr<yap::ParticleCombination> a = yap::ParticleCombination::uniqueSharedPtr(0);
    std::shared_ptr<yap::ParticleCombination> b = yap::ParticleCombination::uniqueSharedPtr(1);
    std::shared_ptr<yap::ParticleCombination> c = yap::ParticleCombination::uniqueSharedPtr(2);
    std::shared_ptr<yap::ParticleCombination> d = yap::ParticleCombination::uniqueSharedPtr(3);

    std::shared_ptr<yap::ParticleCombination> a_b = yap::ParticleCombination::uniqueSharedPtr({a, b});

    std::shared_ptr<yap::ParticleCombination> b2 = yap::ParticleCombination::uniqueSharedPtr(1);

    std::shared_ptr<yap::ParticleCombination> a_b2 = yap::ParticleCombination::uniqueSharedPtr({a, b2});

    std::cout << "b  = " << b << "\nb2 = " << b2 << std::endl;
    std::cout << "a_b  = " << a_b << "\na_b2 = " << a_b2 << std::endl;

    for (auto d : yap::ParticleCombination::particleCombinationSet()) {
        for (yap::ParticleIndex i : d->indices())
            std::cout << static_cast<unsigned>(i) << std::flush;
        std::cout << std::endl;
    }

    // ab->addDaughter(b);

    // std::shared_ptr<yap::ParticleCombination> c_d = std::make_shared<yap::ParticleCombination>();
    // cd->addDaughter(c);
    // cd->addDaughter(d);

    // std::shared_ptr<yap::ParticleCombination> a_d = std::make_shared<yap::ParticleCombination>();
    // ad->addDaughter(a);
    // ad->addDaughter(d);

    // std::shared_ptr<yap::ParticleCombination> c_b = std::make_shared<yap::ParticleCombination>();
    // cb->addDaughter(c);
    // cb->addDaughter(b);

    // std::shared_ptr<yap::ParticleCombination> abc = std::make_shared<yap::ParticleCombination>();
    // abc->addDaughter(ab);
    // abc->addDaughter(c);

    // std::shared_ptr<yap::ParticleCombination> ab_cd = std::make_shared<yap::ParticleCombination>();
    // ab_cd->addDaughter(ab);
    // ab_cd->addDaughter(cd);

    // std::shared_ptr<yap::ParticleCombination> cd_ab = std::make_shared<yap::ParticleCombination>();
    // cd_ab->addDaughter(cd);
    // cd_ab->addDaughter(ab);

    // std::shared_ptr<yap::ParticleCombination> abc_d = std::make_shared<yap::ParticleCombination>();
    // abc_d->addDaughter(abc);
    // abc_d->addDaughter(d);

    // std::shared_ptr<yap::ParticleCombination> ab_ad = std::make_shared<yap::ParticleCombination>();
    // ab_ad->addDaughter(ab);
    // ab_ad->addDaughter(ad);

    // std::cout << "ab_cd = " << std::flush;
    // for (unsigned i = 0; i < ab_cd->indices().size(); ++i)
    //     std::cout << static_cast<unsigned>(ab_cd->indices()[i]) << std::flush;
    // std::cout << std::endl;

    // std::cout << "cd_ab = " << std::flush;
    // for (unsigned i = 0; i < cd_ab->indices().size(); ++i)
    //     std::cout << static_cast<unsigned>(cd_ab->indices()[i]) << std::flush;
    // std::cout << std::endl;

    // std::cout << "ab_ad = " << std::flush;
    // for (unsigned i = 0; i < ab_ad->indices().size(); ++i)
    //     std::cout << static_cast<unsigned>(ab_ad->indices()[i]) << std::flush;
    // std::cout << std::endl;

    // std::cout << (abc_d == ab_cd) << std::endl;

    // std::cout << "ab_ad::consistent() = " << ab_ad->consistent() << std::endl;

    // for (unsigned i = 0; i < D1->daughters().size(); ++i) {
    //     std::cout << i << std::flush;
    //     std::shared_ptr<yap::ParticleCombination> d = D1->daughters()[i].lock();
    //     for (unsigned j = 0; j < d->indices().size(); ++j)
    //         std::cout << "\t" << static_cast<unsigned>(d->indices()[j]) << std::flush;
    //     std::cout << std::endl;
    // }

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
