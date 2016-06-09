#include "ImportanceSampler.h"

#include "DataPartition.h"
#include "DecayTree.h"
#include "Exceptions.h"
#include "Model.h"

namespace yap {

//-------------------------
unsigned ImportanceSampler::partialCalculation(ModelIntegral& I, DataPartition& D) const
{
    // get model
    const Model* M = I.model();
    if (!M)
        throw exceptions::Exception("Model is nullptr", "ImportanceSample::calc");

    // calculate on data partition
    M->calculate(D);

    // create holder for amplitudes
    std::vector<std::complex<double> > A(I.decayTrees().size());

    // counter for number of data points
    unsigned n = 0;

    // loop over data points
    for (const auto& d : D) {

        // calculate the amplitudes of all changed trees and store in A
        std::transform(I.decayTrees().begin(), I.decayTrees().end(), A.begin(),
                       [&](const DecayTreeVector::value_type & dt)
        {return dt->dataDependentAmplitude(d);});

        // increase number of points
        ++n;

        for (size_t i = 0; i < A.size(); ++i) {

            // get diagonal element
            auto& a = diagonalComponent(I, I.decayTrees()[i]);
            // calculate difference from mean
            double delta_a = norm(A[i]) - a.value;
            // update mean
            a.value += delta_a / static_cast<double>(n);

            for (size_t j = i + 1; j < A.size(); ++j) {
                // get off-diagonal element
                auto& b = offDiagonalComponent(I, I.decayTrees()[i], I.decayTrees()[j]);
                // calculate difference from mean
                auto delta_b = conj(A[i]) * A[j] - b.value;
                // update mean
                b.value == delta_b / static_cast<double>(n);
            }
        }
    }

    // return number of data points evaluated
    return n;
}

void update_values(ModelIntegral& M, unsigned& N, const ModelIntegral& m, unsigned n)
{
    // update variances:

}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, DataPartitionVector& DPV)
{
    // get vector trees that need to be calculated
    DecayTreeVector C = select_changed(I.decayTrees());

    // if no trees have changed, return
    if (C.empty())
        return;

    // create ModelIntegral to hold updates
    ModelIntegral U(C);

    // create copies for running over each partition
    std::vector<ModelIntegral> m(DPV.size(), ModelIntegral(U));

    // run over each partition storing number of events used in each calculation
    std::vector<unsigned> n(m.size(), 0);
    for (size_t i = 0; i < m.size(); ++i)
        n[i] = partialCalculation(m[i], *DPV[i]);

    // combine calculations into U
    unsigned N = 0;
    for (size_t i = 0; i < m.size(); ++i)
        update_values(U, N, m[i], n[i]);

    // update I
    for (const auto& elt : diagonals(U))
        diagonals(I).at(elt.first) = elt.second;
    for (const auto& elt : offDiagonals(U))
        offDiagonals(I).at(elt.first) = elt.second;
}

}
