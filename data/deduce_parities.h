#ifndef yap_deduce_parities_h
#define yap_deduce_parities_h

#include <Exceptions.h>
#include <MathUtilities.h>
#include <ParticleTable.h>
#include <QuantumNumbers.h>

/// \return whether pdg code is for a meson
inline bool is_meson(int pdg)
{
    // fundamental particles
    if (abs(pdg) < 100)
        return false;
    // generic mesons
    if (yap::is_odd(abs(pdg) % 10))
        return true;
    // K0L and K0S
    if (abs(pdg) == 310 or abs(pdg) == 130)
        return true;
    return false;
}

/// Deduce quantum number from PDG number;
/// charge & parities assigned with assumption
/// of standard quark quantum numbers.
///
/// Based on scheme given in PDG book,
/// chapter "Monte Carlo Particle Numbering Scheme"
inline int deduce_meson_parity(int pdg)
{
    if (!is_meson(pdg))
        throw yap::exceptions::Exception("pdg is not for a meson", "deduce_meson_parity");
    
    // first check for K0
    if (abs(pdg) == 310 or abs(pdg) == 130)
        return -1;
    
    // break into digits
    std::vector<unsigned> n;
    unsigned p = abs(pdg);
    while (p != 0) {
        n.push_back(p % 10);
        p /= 10;
    }

    // check for exceptions
    // 9___22_:
    if (n.size() > 6 and n[6] == 9 and n[1] == 2 and n[2] == 2)
        throw yap::exceptions::Exception("cannot deduce parity for state with code 9xxx22x", "deduce_meson_parity");

    // check for mistakes in PDG scheme
    // a_0(980) --- should be 9010q11
    if (abs(pdg) == 9000111 or abs(pdg) == 9000211)
        return +1;
    // pi(1800) --- should be 9000q11
    if (abs(pdg) == 9010111 or abs(pdg) == 9010211)
        return -1;
    
    /////////////////////////
    // use encoding of PDG numbers

    // first digit is 2J + 1
    unsigned twoJ = n[0] - 1;

    // orbital angular momenta
    unsigned L = (twoJ == 0) ? 0 : twoJ / 2 - 1;
    if (n.size() > 4) {
        if (twoJ == 0)
            L = n[4];
        else {
            switch (n[4]) {
            case 0:
                L = twoJ / 2 - 1;
                break;
            case 1:
            case 2:
                L = twoJ / 2;
                break;
            case 3:
                L = twoJ / 2 + 1;
                break;
            default:
                if (n.size() < 7 or n[6] != 9)
                    throw yap::exceptions::Exception("n_L > 3 in " + std::to_string(pdg), "deduce_quantum_numbers");
            }
        }
    }
    
    // P = (-1)^(L + 1)
    return yap::pow_negative_one(L + 1);
}

/// Deduce meson parities from PDG numbers
/// \param T ParticleTable to alter
inline void deduce_meson_parities(yap::ParticleTable& T)
{
    for (auto& pdg_t : T) {
        // only correct mesons
        if (!is_meson(pdg_t.first))
            continue;
        try {
            auto P = deduce_meson_parity(pdg_t.first);
            pdg_t.second.quantumNumbers().setP(P);
        } catch (yap::exceptions::Exception) { /* do nothing */ }
    }
}

#endif
