#include <Attributes.h>
#include <ConstantWidthBreitWigner.h>
#include <DataPartition.h>
#include <DataPoint.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <Group.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <MassAxes.h>
#include <MassShapeWithNominalMass.h>
#include <Model.h>
#include <Parameter.h>
#include <Particle.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <PHSP.h>
#include <Sort.h>
#include <SpinAmplitudeCache.h>
#include <WignerD.h>

#include <memory>
#include <random>
#include <string>

//#include <callgrind.h>

int main( int argc, char** argv)
{

    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    yap::Model M(std::make_unique<yap::HelicityFormalism>());

    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : "./data") + "/evt.pdl");
    
    // add some missing parities
    T[421].setQuantumNumbers(yap::QuantumNumbers(0, 0, -1)); // D0
    T[211].setQuantumNumbers(yap::QuantumNumbers(1, 0, -1)); // pi+
    T[113].setQuantumNumbers(yap::QuantumNumbers(0, 2, -1)); // rho0
    T[20213].setQuantumNumbers(yap::QuantumNumbers(1, 2, +1)); // a_1+
    T[9000221].setQuantumNumbers(yap::QuantumNumbers(0, 0, +1)); // sigma

    double radialSize = 1.;

    // initial state particle
    auto D = yap::DecayingParticle::create(T[421], radialSize);

    auto D_mass = T[421].mass();

    // final state particles
    auto piPlus  = yap::FinalStateParticle::create(T[211]);
    auto piMinus = yap::FinalStateParticle::create(T[-211]);

    // Set final-state particles
    M.setFinalState(piPlus, piMinus, piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = yap::DecayingParticle::create(T[9000221], radialSize, std::make_shared<yap::ConstantWidthBreitWigner>(T[9000221]));
    sigma->addStrongDecay(piPlus, piMinus);

    // rho
    auto rho = yap::DecayingParticle::create(T[113], radialSize, std::make_shared<yap::ConstantWidthBreitWigner>(T[113]));
    rho->addStrongDecay(piPlus, piMinus);

    // f_2
    auto f_2 = yap::DecayingParticle::create(T["f_2"], radialSize, std::make_shared<yap::ConstantWidthBreitWigner>(T["f_2"]));
    f_2->addStrongDecay(piPlus, piMinus);
    
    // // omega
    // auto omega = yap::DecayingParticle::create(T[223], radialSize, std::make_shared<yap::ConstantWidthBreitWigner>(T[223]));
    // omega->addStrongDecay(piPlus, piMinus);

    // a_1
    auto a_1 = yap::DecayingParticle::create(T[20213], radialSize, std::make_shared<yap::ConstantWidthBreitWigner>(T[20213]));
    // a_1->addStrongDecay(sigma, piPlus);
    a_1->addStrongDecay(rho,   piPlus);

    // D's channels
    D->addWeakDecay(rho, rho);
    // D->addWeakDecay(omega, omega);
    // D->addWeakDecay(rho, omega);
    D->addWeakDecay(a_1, piMinus);
    D->addWeakDecay(sigma, piPlus, piMinus);
    D->addWeakDecay<2>(f_2, f_2);
    D->addWeakDecay(piPlus, piMinus, piPlus, piMinus);

    // add other background particles
    M.addInitialState(a_1);
    M.addInitialState(rho);

    M.lock();

    // check consistency
    if (!M.consistent()) {
        LOG(INFO) << "inconsistent!";
        return 1;
    }

    // print stuff

    // for (const auto& isp : M.initialStates()) {
    //     FLOG(INFO) << "";
    //     FLOG(INFO) << isp->particleCombinations().size() << " " << *isp << " symmetrizations";
    //     for (auto& pc : isp->particleCombinations())
    //         FLOG(INFO) << *pc;
    //     FLOG(INFO) << "";
    // }
    
    // FLOG(INFO) << "\nFour momenta symmetrizations with " << M.fourMomenta()->nSymmetrizationIndices() << " indices";
    // for (auto& pc_i : M.fourMomenta()->symmetrizationIndices())
    //     FLOG(INFO) << *pc_i.first << ": " << pc_i.second;

    MULTILINE(FLOG(INFO),to_decay_string(*D));
    FLOG(INFO) << "";
    
    MULTILINE(FLOG(INFO),to_string(*M.spinAmplitudeCache()));

    LOG(INFO) << "create dataPoints";

    // create default mass axes
    
    auto A = M.massAxes();
    auto m2r = yap::squared(yap::mass_range(D_mass, A, M.finalStateParticles()));

    // create data set
    yap::DataSet data = M.createDataSet();
    // create random number engine for generation of points
    std::mt19937 g(0);
    // fill data set with 1 point
    std::generate_n(std::back_inserter(data), 1,
                    std::bind(yap::phsp<std::mt19937>, std::cref(M), D_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));

    LOG(INFO) << "Size of DataPoint: " + std::to_string(data[0].bytes()) + " byte (for " + std::to_string(data[0].nDataAccessors()) + " data accessors)";
    LOG(INFO) << "Printing data:";
    for (unsigned d = 0; d < data.size(); ++d) {
        LOG(INFO) << "  DataPoint " << d;
        for (auto& v : M.fourMomenta()->finalStateMomenta(data[d]))
            LOG(INFO) << yap::to_string(v);
    }

    // create data partitions
    unsigned nChains = 1;
    auto parts = yap::DataPartitionWeave::create(data, nChains);

    //CALLGRIND_START_INSTRUMENTATION

    // create uniform random distributions
    std::uniform_real_distribution<double> uniform;
    std::uniform_real_distribution<double> uniform2(0.95, 1.052631579);

    yap::mass_parameter get_mass;
    
    // do several loops over all dataPartitions
    for (unsigned i = 0; i < 10; ++i) {

        // change amplitudes
        if (uniform(g) > 0.5)
            for (auto& isp : M.initialStates())
                for (auto& dt : decay_trees(*isp, yap::is_not_fixed()))
                    if (uniform(g) > 0.5)
                        *dt->freeAmplitude() = uniform2(g) * dt->freeAmplitude()->value();

        // change masses
        if (uniform(g) > 0.5)
            for (auto& d : particles(M, yap::has_a_mass_shape(), yap::has_a_mass()))
                if (uniform(g) > 0.5 and get_mass(*d).variableStatus() != yap::VariableStatus::fixed) {
                    DEBUG("change mass for " << to_string(*d));
                    const_cast<yap::RealParameter&>(get_mass(*d)) = uniform2(g) * get_mass(*d).value();
                }

        DEBUG("===================================================================================================================== ");

        double logA = sum_of_log_intensity(M, parts);
        M.setParameterFlagsToUnchanged();

        LOG(INFO) << "logA = " << logA;

    }

    //CALLGRIND_STOP_INSTRUMENTATION

    /*
        for (auto& a : freeAmps)
            a->setValue(0.5 * a->value());

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);

        // only change some amps
        unsigned i(0);
        for (auto& a : freeAmps) {
            if (i++ % 2 == 0)
                a->setValue(1.);
        }

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);

        // set to zero
        for (auto& a : freeAmps) {
            a->setValue(0.);
        }

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);
    */


    MULTILINE(LOG(INFO), to_string(D->decayTrees()));

    LOG(INFO) << std::endl << "Free amplitudes: ";
    for (const auto& fa : sort(free_amplitudes(M), yap::compare_by<yap::is_fixed>(),
                               yap::by_parent_name<>(), yap::by_l<>()))
        LOG(INFO) << yap::to_string(*fa);

    LOG(INFO) << "";
    LOG(INFO) << "Grouped by l,s";
    for (const auto& fav : group(free_amplitudes(M), yap::by_l<>(), yap::by_s<>())) {
        LOG(INFO) << "" ;
        LOG(INFO) << "L = " << yap::orbital_angular_momentum()(fav[0])
                  << " S = " << yap::spin_angular_momentum()(fav[0])
                  << " :: " << std::to_string(fav.size()) << " amplitudes:";
        for (const auto& fa : sort(fav, yap::compare_by<yap::is_fixed>(), yap::by_parent_name<>()))
            LOG(INFO) << yap::to_string(*fa);
    }
    
    FLOG(INFO) << "alright! \n";

    return 0;
}
