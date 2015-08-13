#include "logging.h"
INITIALIZE_EASYLOGGINGPP

#include "DataSet.h"
#include "DataPoint.h"

#include <TLorentzVector.h>

int main( int argc, char** argv)
{

    std::vector<TLorentzVector> P;
    P.push_back(TLorentzVector(0, 1, 2, 3));
    P.push_back(TLorentzVector(1, 2, 3, 0));
    P.push_back(TLorentzVector(2, 3, 0, 1));
    P.push_back(TLorentzVector(3, 0, 1, 2));

    yap::DataSet DS;
    DS.addDataPoint(yap::DataPoint(P)); // move

    yap::DataPoint D(P);
    DS.addDataPoint(D);         // copy and move
    DS.addDataPoint(std::move(D)); // move

    P.pop_back();
    DS.addDataPoint(yap::DataPoint(P)); // fail

}
