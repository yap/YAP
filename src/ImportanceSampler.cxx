#include "ImportanceSampler.h"

#include "DataPartition.h"
#include "DecayTree.h"
#include "DecayTreeVectorIntegral.h"
#include "Exceptions.h"
#include "IntegralElement.h"
#include "Model.h"
#include "ModelIntegral.h"

#include "logging.h"

#include <future>
#include <vector>

namespace yap {

//-------------------------
unsigned ImportanceSampler::partially_calculate_subIntegral(DecayTreeVectorIntegral& I, DataPartition& D)
{
    // create holder for amplitudes
    std::vector<std::complex<double> > A(I.decayTrees().size());

    // create holder for values
    std::vector<DiagonalIntegralMap::mapped_type*> diags;
    std::vector<std::vector<OffDiagonalIntegralMap::mapped_type*> > off_diags(I.decayTrees().size() - 1);
    // fill holders
    diags.reserve(I.decayTrees().size());
    for (size_t i = 0; i < I.decayTrees().size(); ++i) {
        diags.push_back(&diagonalComponent(I, I.decayTrees()[i]));
        if (!off_diags.empty()) {
            off_diags[i].reserve(I.decayTrees().size() - i - 1);
            for (size_t j = i + 1; j < I.decayTrees().size(); ++j) {
                off_diags[i].push_back(&offDiagonalComponent(I, I.decayTrees()[i], I.decayTrees()[j]));
            }
        }
    }

    // counter for number of data points
    unsigned n = 0;

    // loop over data points
    for (const auto& d : D) {

        // calculate the amplitudes of all changed trees and store in A
        std::transform(I.decayTrees().begin(), I.decayTrees().end(), A.begin(),
                       [&d](const DecayTreeVector::value_type & dt)
        {return dt->dataDependentAmplitude(d);});

        // increase number of points
        ++n;

        for (size_t i = 0; i < A.size(); ++i) {

            // calculate difference from mean
            double delta_diag = norm(A[i]) - diags[i]->value;
            // update mean
            diags[i]->value += delta_diag / static_cast<double>(n);

            for (size_t j = i + 1; j < A.size(); ++j) {
                // calculate difference from mean
                auto delta_offdiag = conj(A[i]) * A[j] - off_diags[i][j - i - 1]->value;
                // update mean
                off_diags[i][j - i - 1]->value += delta_offdiag / static_cast<double>(n);
            }
        }
    }

    // return number of data points evaluated
    return n;
}

//-------------------------
unsigned ImportanceSampler::partially_calculate(integral_sub_map& J, DataPartition& D)
{
    // get model
    const Model* M = J.begin()->second.model();
    if (!M)
        throw exceptions::Exception("Model is nullptr", "ImportanceSampler::partialCalculate");

    // calculate on data partition
    M->calculate(D);

    unsigned n = 0;
    for (auto& j : J)
        n = ImportanceSampler::partially_calculate_subIntegral(j.second, D);

    return n;
}

//-------------------------
ImportanceSampler::integral_sub_map ImportanceSampler::select_changed_integrals(ModelIntegral& I)
{
    ImportanceSampler::integral_sub_map J;

    for (auto& b_I : integrals(I)) {

        // get vector trees that need to be calculated
        auto C = select_changed(b_I.second.decayTrees());

        if (!C.empty())
            J.emplace(&b_I.second, DecayTreeVectorIntegral(C));

    }

    return J;
}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, DataPartition& D)
{
    // get DecayTreeVectorIntegral's for DecayTree's that need to be calculated
    auto J = select_changed_integrals(I);

    // if no trees have changed, return
    if (J.empty())
        return;

    // calculate it
    partially_calculate(J, D);

    // copy calculated values
    for (auto& p_i : J) {
        const auto& dtv = p_i.second.decayTrees();
        // set results
        for (size_t i = 0; i < dtv.size(); ++i) {
            // diagonal components
            diagonalComponent(*p_i.first, dtv[i]) = diagonalComponent(p_i.second, dtv[i]);
            // diagonal components
            for (size_t j = i + 1; j < dtv.size(); ++j)
                offDiagonalComponent(*p_i.first, dtv[i], dtv[j]) = offDiagonalComponent(p_i.second, dtv[i], dtv[j]);
        }
    }
}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, DataPartitionVector& DPV)
{
    // if only one data partition, don't thread:
    if (DPV.size() == 1) {
        calculate(I, *DPV[0]);
        return;
    }

    // get DecayTreeVectorIntegral's for DecayTree's that need to be calculated
    auto J = select_changed_integrals(I);

    // if no trees have changed, return
    if (J.empty())
        return;

    // create copies for running over each partition
    std::vector<integral_sub_map> m(DPV.size(), J);

    // run over each partition storing number of events used in each calculation
    std::vector<std::future<unsigned> > n;
    n.reserve(m.size());
    // create thread for each partial calculation
    for (size_t i = 0; i < m.size(); ++i)
        n.push_back(std::async(std::launch::async, &ImportanceSampler::partially_calculate,
                               std::ref(m[i]), std::ref(*DPV[i])));

    // calculate data fractions:
    // also waits for threads to finish calculating
    std::vector<double> f;
    f.reserve(n.size());
    std::transform(n.begin(), n.end(), std::back_inserter(f), std::mem_fn(&std::future<unsigned>::get));
    double N = std::accumulate(f.begin(), f.end(), 0.);
    std::transform(f.begin(), f.end(), f.begin(), std::bind(std::divides<double>(), std::placeholders::_1, N));

    // reset current state for newly calculated results
    for (auto& p_i : J) {
        const auto& dtv = p_i.second.decayTrees();
        for (size_t i = 0; i < dtv.size(); ++i) {
            // reset diagonal component
            diagonalComponent(*p_i.first, dtv[i]).value = 0;
            for (size_t j = i + 1; j < dtv.size(); ++j) {
                // reset off-diagonal component
                offDiagonalComponent(*p_i.first, dtv[i], dtv[j]).value = 0;
            }
        }
    }

    for (size_t k = 0; k < m.size(); ++k) {
        for (auto& p_i : m[k]) {
            const auto& dtv = p_i.second.decayTrees();
            for (size_t i = 0; i < dtv.size(); ++i) {
                // combine diagonal components
                diagonalComponent(*p_i.first, dtv[i]).value += f[k] * diagonalComponent(p_i.second, dtv[i]).value;
                // off-diagonal elements
                for (size_t j = i + 1; j < dtv.size(); ++j) {
                    offDiagonalComponent(*p_i.first, dtv[i], dtv[j]).value += f[k] * offDiagonalComponent(p_i.second, dtv[i], dtv[j]).value;
                }
            }
        }
    }
}

}
