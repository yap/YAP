#include "PDL.h"

#include "ParticleFactory.h"

#include <fstream>
#include <sstream>
#include <string>

namespace yap {

//-------------------------
const ParticleTableEntry PDLIterator::operator*() const
{
    // create a stream with the entry string
    std::istringstream iss(Value_);

    // temporary variables
    std::string discard, particleName;
    int    stdhepid;
    double mass;
    double particleWidth;
    double particleMaxWidth; // not needed
    int    threeTimesCharge;
    int    twoTimesSpin;
    double lifetime; // dimension [c/mm]
    int    lundkc; // not needed

    // read the values in
    iss >> discard; // "add"
    iss >> discard; // "p"
    iss >> discard; // "Particle"
    iss >> particleName;
    iss >> stdhepid;
    iss >> mass;
    iss >> particleWidth;
    iss >> particleMaxWidth;
    iss >> threeTimesCharge;
    iss >> twoTimesSpin;
    iss >> lifetime;
    iss >> lundkc;

    return ParticleTableEntry(stdhepid, particleName,
                              QuantumNumbers(twoTimesSpin, std::round(1. * threeTimesCharge / 3)),
                              mass, {particleWidth});
}

//-------------------------
PDLIterator& PDLIterator::operator++()
{
    // ignore comments and 'set' entries
    std::string command;
    do {
        // check if iterator reaches EOF
        if (InputStream_ && !getline(*InputStream_, Value_)) {
            *this = end();
        } else {
            std::istringstream iss(Value_);
            iss >> command;
            // check if the 'end' keyword has been reached
            if (command == "end")
                *this = end();
        }
    } while (Value_[0] == '*' || command == "set");

    return *this;
}

//-------------------------
PDLIterator PDLIterator::operator++(int)
{
    PDLIterator prev(*this);
    // read line in
    ++(*this);
    return prev;
}

//-------------------------
ParticleFactory read_pdl_file(const std::string& filename)
{
    std::ifstream inFile(filename, std::ios::in);
    ParticleFactory particleFactory;
    std::copy(PDLIterator(inFile), PDLIterator::end(),
              inserter(particleFactory));
    return particleFactory;
}

}
