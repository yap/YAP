#include "ImportanceSampler.h"

#include "DataPartition.h"
#include "DecayTree.h"
#include "Exceptions.h"
#include "Model.h"

#include "logging.h"

#include <future>
#include <vector>

namespace yap {

//-------------------------
unsigned ImportanceSampler::partialCalculation(ModelIntegral& I, DataPartition& D)
{
    // get model
    const Model* M = I.model();
    if (!M)
        throw exceptions::Exception("Model is nullptr", "ImportanceSampler::partialCalculate");

    // calculate on data partition
    M->calculate(D);

    // create holder for amplitudes
    std::vector<std::complex<double> > A(I.decayTrees().size());

    // create holder for values
    std::vector<ModelIntegral::DiagonalMap::mapped_type*> diags;
    std::vector<std::vector<ModelIntegral::OffDiagonalMap::mapped_type*> > off_diags(I.decayTrees().size() - 1);
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
                off_diags[i][j - i - 1]->value = delta_offdiag / static_cast<double>(n);
            }
        }
    }

    // return number of data points evaluated
    return n;
}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, DataPartition& D)
{
    // get vector trees that need to be calculated
    DecayTreeVector C = select_changed(I.decayTrees());

    // if no trees have changed, return
    if (C.empty())
        return;

    // create ModelIntegral to hold updates
    ModelIntegral U(C);
    // calculate it
    partialCalculation(U, D);

    // set results
    for (size_t i = 0; i < C.size(); ++i) {
        // diagonal components
        diagonalComponent(I, C[i]) = diagonalComponent(U, C[i]);
        // diagonal components
        for (size_t j = i + 1; j < C.size(); ++j)
            offDiagonalComponent(I, C[i], C[j]) = offDiagonalComponent(U, C[i], C[j]);
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

    // get vector trees that need to be calculated
    DecayTreeVector C = select_changed(I.decayTrees());

    // if no trees have changed, return
    if (C.empty())
        return;

    // create copies for running over each partition
    std::vector<ModelIntegral> m(DPV.size(), ModelIntegral(C));

    // run over each partition storing number of events used in each calculation
    std::vector<std::future<unsigned> > n;
    n.reserve(m.size());
    // create thread for each partial calculation
    for (size_t i = 0; i < m.size(); ++i)
        n.push_back(std::async(std::launch::async, &ImportanceSampler::partialCalculation,
                               std::ref(m[i]), std::ref(*DPV[i])));

    // calculate data fractions:
    // also waits for threads to finish calculating
    double N = std::accumulate(n.begin(), n.end(), 0., [](double a, std::future<unsigned>& F) {return a + F.get();});
    std::vector<double> f(n.size());
    std::transform(n.begin(), n.end(), f.begin(), [&](std::future<unsigned>& nn) {return nn.get() / N;});

    // combine results
    for (size_t i = 0; i < C.size(); ++i) {
        // reset diagonal component
        diagonalComponent(I, C[i]).value = 0;
        // combine diagonal components
        for (size_t k = 0; k < m.size(); ++k)
            diagonalComponent(I, C[i]).value += f[k] * diagonalComponent(m[k], C[i]).value;
        // off-diagonal elements
        for (size_t j = i + 1; j < C.size(); ++j) {
            // reset off-diagonal component
            offDiagonalComponent(I, C[i], C[j]).value = 0;
            for (size_t k = 0; k < m.size(); ++k)
                offDiagonalComponent(I, C[i], C[j]).value += f[k] * offDiagonalComponent(m[k], C[i], C[j]).value;
        }
    }
}

}
