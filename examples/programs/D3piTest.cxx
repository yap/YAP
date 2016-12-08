#include <Attributes.h>
#include <BreitWigner.h>
#include <container_utils.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <DecayTreeVectorIntegral.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <PoleMass.h>
#include <PHSP.h>
#include <RelativisticBreitWigner.h>
#include <Sort.h>
#include <ZemachFormalism.h>

#include <memory>
#include <random>
#include <vector>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    LOG(INFO) << "Start";

    // yap::Model M(std::make_unique<yap::ZemachFormalism>());
    yap::Model M(std::make_unique<yap::HelicityFormalism>());

    LOG(INFO) << "Model created";

    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    LOG(INFO) << "table loaded";

    // add some missing parities
    T["D+"].setQuantumNumbers(yap::QuantumNumbers(1, 0, -1));
    T[211].setQuantumNumbers(yap::QuantumNumbers(1, 0, -1));
    T[113].setQuantumNumbers(yap::QuantumNumbers(0, 2, -1));
    T[225].setQuantumNumbers(yap::QuantumNumbers(0, 4, +1));
    T["f_0"].setQuantumNumbers(yap::QuantumNumbers(0, 0, +1));
    T["f_0(500)"].setQuantumNumbers(T["f_0"].quantumNumbers());
    T["f_0(1500)"].setQuantumNumbers(T["f_0"].quantumNumbers());
        
    // initial state particle
    auto D = yap::DecayingParticle::create(T["D+"], radialSize);

    auto D_mass = T["D+"].mass();

    LOG(INFO) << "D created";

    // final state particles
    auto piPlus  = yap::FinalStateParticle::create(T[211]);
    auto piMinus = yap::FinalStateParticle::create(T[-211]);

    LOG(INFO) << "fsp's created";

    // set final state
    M.setFinalState(piPlus, piMinus, piPlus);

    LOG(INFO) << "final state set";

    // rho
    auto rho = yap::DecayingParticle::create(T[113], radialSize, std::make_shared<yap::RelativisticBreitWigner>(T[113]));
    rho->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(rho, piPlus);

    // f_2(1270)
    auto f_2 = yap::DecayingParticle::create(T[225], radialSize, std::make_shared<yap::BreitWigner>(T[225]));
    f_2->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_2, piPlus);

    // f_0(980)
    auto f_0_980_flatte = std::make_shared<yap::Flatte>(0.965);
    f_0_980_flatte->add(yap::FlatteChannel(0.406, *piPlus, *piMinus));
    f_0_980_flatte->add(yap::FlatteChannel(0.406 * 2, T[321], T[-321])); // K+K-
    auto f_0_980 = yap::DecayingParticle::create("f_0_980", T["f_0"].quantumNumbers(), radialSize, f_0_980_flatte);
    f_0_980->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_0_980, piPlus);

    // f_0(1370)
    auto f_0_1370 = yap::DecayingParticle::create("f_0_1370", T["f_0"].quantumNumbers(), radialSize, std::make_unique<yap::BreitWigner>(1.350, 0.265));
    f_0_1370->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_0_1370, piPlus);

    // f_0(1500)
    auto f_0_1500 = yap::DecayingParticle::create(T["f_0(1500)"], radialSize, std::make_unique<yap::BreitWigner>(T["f_0(1500)"]));
    f_0_1500->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_0_1500, piPlus);

    // sigma a.k.a. f_0(500)
    auto sigma = yap::DecayingParticle::create(T["f_0(500)"], radialSize, std::make_unique<yap::PoleMass>(std::complex<double>(0.470, -0.220)));
    sigma->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(sigma, piPlus);

    FLOG(INFO) << "number of free amplitudes = " << free_amplitudes(M).size();

    // Add channels to D
    *free_amplitude(*D, yap::to(rho))      = std::polar(1., 0.);
    *free_amplitude(*D, yap::to(f_0_980))  = std::polar(1.4, yap::rad(12.));
    *free_amplitude(*D, yap::to(f_2))      = std::polar(2.1, yap::rad(-123.));
    *free_amplitude(*D, yap::to(f_0_1370)) = std::polar(1.3, yap::rad(-21.));
    *free_amplitude(*D, yap::to(f_0_1500)) = std::polar(1.1, yap::rad(-44.));
    *free_amplitude(*D, yap::to(sigma))    = std::polar(3.7, yap::rad(-3.));
    // D->addWeakDecay(piPlus, piMinus, piPlus);

    // check consistency
    if (M.consistent())
        LOG(INFO) << "consistent!";
    else
        LOG(INFO) << "inconsistent!";

    LOG(INFO) << "before locking:";
    for (const auto& dt : rho->decayTrees())
        LOG(INFO) << to_string(*dt);
    
    M.lock();

    LOG(INFO) << "after locking:";
    for (const auto& dt : rho->decayTrees())
        LOG(INFO) << to_string(*dt);

    // print stuff
 
    FLOG(INFO) << "";
    FLOG(INFO) << D->particleCombinations().size() << " D symmetrizations";

    FLOG(INFO) << "";
    FLOG(INFO) << "Four momenta symmetrizations with " << M.fourMomenta()->nSymmetrizationIndices() << " indices";

    // FLOG(INFO) << "\nHelicity angle symmetrizations with " << M.helicityAngles()->maxSymmetrizationIndex() + 1 << " indices \n";

    MULTILINE(FLOG(INFO),to_decay_string(*D));
    FLOG(INFO) << "";

    MULTILINE(FLOG(INFO),to_string(*M.spinAmplitudeCache()));

    // get default Dalitz axes
    auto A = M.massAxes();
    auto m2r = yap::squared(yap::mass_range(D_mass, A, M.finalStateParticles()));

    // generate points randomly in phase space of model
    std::mt19937 g(0);

    // create data set
    auto data = M.createDataSet();

    // generate 10,000 phase-space-distributed data points
    std::generate_n(std::back_inserter(data), 10000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(M), D_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));

    LOG(INFO) << data.size() << " data points of " << data[0].bytes() << " bytes each = " << data.bytes() * 1.e-6 << " MB";

    M.calculate(data);

    yap::ModelIntegral MI(M);
    yap::ImportanceSampler::calculate(MI, data);

    for (const auto& mci : MI.integrals()) {

        MULTILINE(LOG(INFO), to_string(mci.Integral.decayTrees()));

        auto A_DT = amplitude(mci.Integral.decayTrees(), data[0]);
        LOG(INFO) << "A_DT = " << A_DT;
        LOG(INFO) << "|A_DT|^2 = " << norm(A_DT);

        LOG(INFO) << "integral = " << to_string(integral(mci.Integral));
        auto ff = fit_fractions(mci.Integral);
        for (size_t i = 0; i < ff.size(); ++i)
            LOG(INFO) << "fit fraction " << to_string(100. * ff[i]) << "% for " << to_string(*mci.Integral.decayTrees()[i]->freeAmplitude());
        LOG(INFO) << "sum of fit fractions = " << to_string(std::accumulate(ff.begin(), ff.end(), yap::RealIntegralElement()));

        LOG(INFO) << "cached integral components:";
        auto I_cached = cached_integrals(mci.Integral);
        for (const auto& row : I_cached)
            LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                         [](std::string & s, const yap::ComplexIntegralElement & c)
        { return s += "\t" + to_string(c);}).erase(0, 1);

        LOG(INFO) << "integral components:";
        auto I = integrals(mci.Integral);
        for (const auto& row : I)
            LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                         [](std::string & s, const yap::ComplexIntegralElement & c)
                                         { return s += "\t" + to_string(c);}).erase(0, 1);
    }

    LOG(INFO) << std::endl << "Fixed amplitudes: ";
    for (const auto& fa : free_amplitudes(M, yap::is_fixed()))
        LOG(INFO) << yap::to_string(*fa);

    LOG(INFO) << std::endl << "Free amplitudes: ";
    for (const auto& fa : free_amplitudes(M, yap::is_not_fixed()))
        LOG(INFO) << yap::to_string(*fa);

    for (unsigned l = 0; l <= 2; ++l) {
        LOG(INFO) << std::endl << "Amplitudes with l = " << l << ": ";
        for (const auto& fa : free_amplitudes(M, yap::l_equals(l)))
            LOG(INFO) << yap::to_string(*fa);
    }

    LOG(INFO) << std::endl << "Free amplitudes (sorted by fixed/notfixed, parent_name, l):";
    for (const auto& fa : sort(free_amplitudes(M), yap::compare_by<yap::is_fixed>(),
                               yap::by_parent_name<>(), yap::by_l<>()))
        LOG(INFO) << yap::to_string(*fa);
    
    LOG(INFO) << "alright!";

    return 0;
}
