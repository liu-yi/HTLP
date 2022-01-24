#include "HTLP.hpp"

HTLP::HTLP(const long modulus_len, const long T, const long kappa, bool cheeting_mode) : modulus_len_(modulus_len), T_(T), kappa_(kappa), cheeting_mode_(cheeting_mode)
{
    rsa_ = GenerateRSAModulus(modulus_len);
    n_ = rsa_.p * rsa_.q;
    n_square_ = n_ * n_;
    lambda_ = (rsa_.p - 1) * (rsa_.q - 1) / 2;
    trapdoor_ = PowerMod(NTL::ZZ(2), T_, lambda_);
    g_ = GenerateJacobiOne();
    h_ = NTL::PowerMod(g_, trapdoor_, n_);
    prime_len_ = ceil((2.0 * kappa_ * log(2) - log(2 * kappa_ * log(2) - 1.1)) / 8); // For x-th prime y, we have y < x / (log(x) - 1.1)
}

HTLP::HTLP(const long modulus_len, const long T, const long kappa) : modulus_len_(modulus_len), T_(T), kappa_(kappa), cheeting_mode_(true)
{
    rsa_ = GenerateRSAModulus(modulus_len);
    n_ = rsa_.p * rsa_.q;
    n_square_ = n_ * n_;
    lambda_ = (rsa_.p - 1) * (rsa_.q - 1) / 2;
    trapdoor_ = PowerMod(NTL::ZZ(2), T_, lambda_);
    g_ = GenerateJacobiOne();
    h_ = NTL::PowerMod(g_, trapdoor_, n_);
}

HTLP::HTLP(const NTL::ZZ &n, const NTL::ZZ &g, const NTL::ZZ &h, const long T, const long kappa) : modulus_len_(NumBits(n)), g_(g), h_(h), T_(T), kappa_(kappa), cheeting_mode_(false)
{
    n_square_ = n_ * n_;
}

RSA HTLP::GenerateRSAModulus(const long modulus_len)
{
    NTL::ZZ p, q;
    while (true)
    {
        p = 2 * NTL::GenGermainPrime_ZZ(modulus_len / 2 - 1) + 1;
        q = 2 * NTL::GenGermainPrime_ZZ(modulus_len / 2 - 1) + 1;
        n_ = p * q;
        if (NumBits(n_) == modulus_len)
        {
            break;
        }
    }
    rsa_.p = p;
    rsa_.q = q;
    return rsa_;
}

NTL::ZZ HTLP::GenerateJacobiOne()
{
    NTL::ZZ a;
    while (true)
    {
        a = RandomBnd(n_);
        if (Jacobi(a, n_) == 1)
        {
            return a;
        }
    }
}

NTL::ZZ HTLP::HashToElement(const std::string str)
{
    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_CTX sha256;
    SHA256_Init(&sha256);
    SHA256_Update(&sha256, str.c_str(), str.size());
    SHA256_Final(hash, &sha256);
    std::stringstream ss;
    NTL::ZZ hash_value = NTL::ZZFromBytes(hash, kappa_ / 8);
    return hash_value;
}

NTL::ZZ HTLP::HashToPrime(const NTL::ZZ &g, const NTL::ZZ &h)
{
    assert(g >= 0 && g < n_);
    assert(h >= 0 && h < n_);

    uint64_t j = 0;
    unsigned char hash[SHA512_DIGEST_LENGTH];
    while (true)
    {
        SHA512_CTX sha512;
        SHA512_Init(&sha512);
        std::string s = "prime";
        SHA512_Update(&sha512, &j, 8);
        SHA512_Update(&sha512, s.c_str(), s.size());
        unsigned char g_bytes[NumBytes(g)];
        unsigned char h_bytes[NumBytes(h)];
        BytesFromZZ(g_bytes, g, NumBytes(g));
        BytesFromZZ(h_bytes, h, NumBytes(g));
        SHA512_Update(&sha512, g_bytes, NumBytes(g));
        SHA512_Update(&sha512, h_bytes, NumBytes(h));
        SHA512_Final(hash, &sha512);
        NTL::ZZ n = NTL::ZZFromBytes(hash, prime_len_);
        if (ProbPrime(n))
        {
            return n;
        }
        j++;
    }
}

NTL::ZZ HTLP::GenerateProof(const long k, const long gamma, std::vector<NTL::ZZ> &C, NTL::ZZ &l)
{
    long k1 = k >> 2;
    long k0 = k - k1;
    NTL::ZZ x(1);
    long k_exp = 1 << k;
    long k0_exp = 1 << k0;
    long k1_exp = 1 << k1;
    auto GetBlock = ([&](long i) -> long
                     {
        NTL::ZZ p = PowerMod(NTL::ZZ(2), T_ - k * (i + 1), l);
        return trunc_long(k_exp * p / l, sizeof(long)); });

    for (long j = gamma - 1; j >= 0; j--)
    {
        x = PowerMod(x, k_exp, n_);
        std::vector<NTL::ZZ> y(k_exp, NTL::ZZ(1));
        long bound = T_ % (k * gamma) == 0 ? T_ / (k * gamma) : T_ / (k * gamma) + 1;
        for (long i = 0; i < bound; i++)
        {
            if (T_ - k * (i * gamma + j + 1) < 0)
            {
                continue;
            }
            long b = GetBlock(i * gamma + j);
            y[b] = y[b] * C[i] % n_;
        }

        for (long b1 = 0; b1 < NTL::power_long(2, k1); b1++)
        {
            NTL::ZZ z(1);
            for (long b0 = 0; b0 < NTL::power_long(2, k0); b0++)
            {
                z = z * y[b1 * NTL::power_long(2, k0) + b0] % n_;
            }
            x = x * PowerMod(PowerMod(z, b1, n_), NTL::power_long(2, k0), n_) % n_;
        }

        for (long b0 = 0; b0 < NTL::power_long(2, k0); b0++)
        {
            NTL::ZZ z(1);
            for (long b1 = 0; b1 < NTL::power_long(2, k1); b1++)
            {
                z = z * y[b1 * NTL::power_long(2, k0) + b0] % n_;
            }
            x = x * PowerMod(z, b0, n_) % n_;
        }
    }
    return x;
}