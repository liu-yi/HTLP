#ifndef AHTLP_HPP_
#define AHTLP_HPP_

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

class AHTLP : public HTLP
{

public:
    AHTLP(const long modulus_len, const long T, const long kappa);
    AHTLP(const NTL::ZZ &n, const NTL::ZZ &g, const NTL::ZZ &h, const long T, const long kappa);
    AHTLP(const long modulus_len, const long T, const long kappa, bool cheeting_mode);

    APuzzle GeneratePuzzle(const NTL::ZZ &s);
    APuzzle GeneratePuzzle(const NTL::ZZ &s, const NTL::ZZ &r);
    NTL::ZZ SolvePuzzle(const APuzzle &Z);
    NTL::ZZ QuickSolvePuzzle(const APuzzle &Z);

    std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> GenerateAValidProof(const APuzzle &Z, const NTL::ZZ &s, const NTL::ZZ &r);
    bool VerifyAValidProof(const APuzzle &Z, const std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> &proof);

    std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> SolvePuzzleWithProof(const long k, const long gamma, const APuzzle &Z);
    std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> QuickSolvePuzzleWithProof(const APuzzle &Z);

    int VerifyProofOfSol(const APuzzle Z, const std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> &proof);
};

#endif