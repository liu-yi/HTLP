#ifndef MHTLP_HPP_
#define MHTLP_HPP_

#include <NTL/ZZ.h>
#include <assert.h>
#include <openssl/sha.h>
#include <vector>
#include <sstream>

#include "Puzzle.hpp"
#include "HTLP.hpp"

#ifndef RSA_
#define RSA_
typedef struct RSA
{
    NTL::ZZ p;
    NTL::ZZ q;
} RSA;
#endif

class MHTLP : public HTLP
{

private:
    NTL::ZZ chi_;

public:
    MHTLP(const long modulus_len, const long T, const long kappa);
    MHTLP(const NTL::ZZ &n, const NTL::ZZ &g, const NTL::ZZ &h, const NTL::ZZ &chi, const long T, const long kappa);
    MHTLP(const long modulus_len, const long T, const long kappa, bool cheeting_mode);

    NTL::ZZ chi()
    {
        return chi_;
    }

    MPuzzle GeneratePuzzle(const NTL::ZZ &s);
    MPuzzle GeneratePuzzle(const NTL::ZZ &s, const NTL::ZZ &r, const NTL::ZZ &r_prime);
    NTL::ZZ SolvePuzzle(const MPuzzle &Z);
    NTL::ZZ QuickSolvePuzzle(const MPuzzle &Z);

    std::tuple<std::vector<NTL::ZZ>, std::vector<NTL::ZZ>, std::vector<NTL::ZZ>, std::vector<NTL::ZZ>> GenerateMValidProof(const MPuzzle &Z, const NTL::ZZ &s, const NTL::ZZ &r, const NTL::ZZ &r_prime);
    bool VerifyMValidProof(const MPuzzle &Z, const std::tuple<std::vector<NTL::ZZ>, std::vector<NTL::ZZ>, std::vector<NTL::ZZ>, std::vector<NTL::ZZ>> &proof);

    std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> SolvePuzzleWithProof(const long k, const long gamma, const MPuzzle &Z);
    std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> QuickSolvePuzzleWithProof(const MPuzzle &Z);

    int VerifyProofOfSol(const MPuzzle Z, const std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> &proof);
};

#endif