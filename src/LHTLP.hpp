#ifndef LHTLP_HPP_
#define LHTLP_HPP_

#include <NTL/ZZ.h>
#include <assert.h>
#include <openssl/sha.h>
#include <vector>
#include <sstream>

#include "HTLP.hpp"
#include "Puzzle.hpp"

#ifndef RSA_
#define RSA_
typedef struct RSA
{
    NTL::ZZ p;
    NTL::ZZ q;
} RSA;
#endif

class LHTLP : public HTLP
{

public:
    LHTLP(const long modulus_len, const long T, const long kappa);
    LHTLP(const NTL::ZZ &n, const NTL::ZZ &g, const NTL::ZZ &h, const long T, const long kappa);
    LHTLP(const long modulus_len, const long T, const long kappa, bool cheeting_mode);

    LPuzzle GeneratePuzzle(const NTL::ZZ &s);
    LPuzzle GeneratePuzzle(const NTL::ZZ &s, const NTL::ZZ &r);
    NTL::ZZ SolvePuzzle(const LPuzzle &Z);
    NTL::ZZ QuickSolvePuzzle(const LPuzzle &Z);

    std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> GenerateLValidProof(const LPuzzle &Z, const NTL::ZZ &s, const NTL::ZZ &r);
    bool VerifyLValidProof(const LPuzzle &Z, const std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> &proof);

    std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> SolvePuzzleWithProof(const long k, const long gamma, const LPuzzle &Z);
    std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> QuickSolvePuzzleWithProof(const LPuzzle &Z);

    int VerifyProofOfSol(const LPuzzle Z, const std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> &proof);
};

#endif