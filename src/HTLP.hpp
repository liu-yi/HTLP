#ifndef HTLP_HPP_
#define HTLP_HPP_

#include <NTL/ZZ.h>
#include <assert.h>
#include <openssl/sha.h>
#include <vector>
#include <sstream>

#include "Puzzle.hpp"

#ifndef RSA_
#define RSA_
typedef struct RSA
{
    NTL::ZZ p;
    NTL::ZZ q;
} RSA;
#endif

class HTLP
{
protected:
    const bool cheeting_mode_;
    RSA rsa_;
    NTL::ZZ n_;
    NTL::ZZ n_square_;
    NTL::ZZ g_;
    NTL::ZZ h_;
    NTL::ZZ lambda_;
    const long T_;
    const long kappa_;
    long prime_len_;
    const long modulus_len_;
    NTL::ZZ trapdoor_;
    RSA GenerateRSAModulus(const long modulus_len);

public:
    HTLP(const long modulus_len, const long T, const long kappa);
    HTLP(const NTL::ZZ &n, const NTL::ZZ &g, const NTL::ZZ &h, const long T, const long kappa);
    HTLP(const long modulus_len, const long T, const long kappa, bool cheeting_mode);

    NTL::ZZ HashToElement(const std::string str);
    NTL::ZZ HashToPrime(const NTL::ZZ &g, const NTL::ZZ &h);

    NTL::ZZ GenerateProof(const long k, const long gamma, std::vector<NTL::ZZ> &C, NTL::ZZ &l);

    NTL::ZZ GenerateJacobiOne();
    NTL::ZZ GenerateRandomExponent()
    {
        return RandomBnd(n_ / 2);
    }
    NTL::ZZ GenerateRandomElement()
    {
        return RandomBnd(n_);
    }

    NTL::ZZ n()
    {
        return n_;
    }
    NTL::ZZ n_square()
    {
        return n_square_;
    }
    NTL::ZZ g()
    {
        return g_;
    }
    NTL::ZZ h()
    {
        return h_;
    }
    long T()
    {
        return T_;
    }

};

#endif