#include "PDL.h"

#include "ParticleTable.h"
#include "set_parities.h"

#include <cmath>
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
                              QuantumNumbers(std::round(1. * threeTimesCharge / 3), twoTimesSpin),
                              mass, {particleWidth});
}

//-------------------------
PDLIterator& PDLIterator::operator++()
{
    // ignore lines not starting with "add"
    std::string command;
    do {

        // check if iterator reaches EOF
        if (InputStream_ and !getline(*InputStream_, Value_))
            return *this = end();

        std::istringstream iss(Value_);
        iss >> command;

        // check if the 'end' keyword has been reached
        if (command == "end")
            return *this = end();

    } while (command != "add");

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
ParticleTable read_pdl_file(const std::string& filename)
{
    std::ifstream inFile(filename, std::ios::in);
    ParticleTable particleFactory;
    std::copy(PDLIterator(inFile), PDLIterator::end(),
              inserter(particleFactory));
    set_parities(particleFactory);
    return particleFactory;
}

}
