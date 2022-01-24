#include "LHTLP.hpp"

#define IS_CHEET_MODE assert(cheeting_mode_);

LHTLP::LHTLP(const long modulus_len, const long T, const long kappa) : HTLP(modulus_len, T, kappa)
{
}
LHTLP::LHTLP(const NTL::ZZ &n, const NTL::ZZ &g, const NTL::ZZ &h, const long T, const long kappa) : HTLP(n, g, h, T, kappa)
{
}
LHTLP::LHTLP(const long modulus_len, const long T, const long kappa, bool cheeting_mode) : HTLP(modulus_len, T, kappa, cheeting_mode)
{
}

std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> LHTLP::GenerateLValidProof(const LPuzzle &Z, const NTL::ZZ &s, const NTL::ZZ &r)
{
    NTL::ZZ x = NTL::RandomLen_ZZ(modulus_len_ - 1 + 2 * kappa_);
    NTL::ZZ t = GenerateRandomElement();
    NTL::ZZ a = PowerMod(g_, x, n_);
    NTL::ZZ b = PowerMod(h_, x * n_, n_square_) * (1 + t * n_) % n_square_;

    std::stringstream ss;
    ss << n_ << g_ << h_ << Z.u << Z.v << a << b;
    NTL::ZZ e = HashToElement(ss.str());

    NTL::ZZ alpha = r * e + x;
    NTL::ZZ beta = s * e + t % n_;

    return {a, b, alpha, beta};
}

bool LHTLP::VerifyLValidProof(const LPuzzle &Z, const std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> &proof)
{
    const NTL::ZZ &a = std::get<0>(proof);
    const NTL::ZZ &b = std::get<1>(proof);
    const NTL::ZZ &alpha = std::get<2>(proof);
    const NTL::ZZ &beta = std::get<3>(proof);
    std::stringstream ss;
    ss << n_ << g_ << h_ << Z.u << Z.v << a << b;
    NTL::ZZ e = HashToElement(ss.str());
    if (PowerMod(Z.u, e, n_) * a % n_ != PowerMod(g_, alpha, n_))
    {
        return false;
    }
    if (PowerMod(h_, alpha * n_, n_square_) * (1 + beta * n_) % n_square_ != PowerMod(Z.v, e, n_square_) * b % n_square_)
    {
        return false;
    }
    return true;
}

LPuzzle LHTLP::GeneratePuzzle(const NTL::ZZ &s)
{
    assert(s < n_ && s >= 0);
    LPuzzle Z(n_);
    NTL::ZZ r = GenerateRandomExponent();
    Z.u = PowerMod(g_, r, n_);
    Z.v = PowerMod(h_, r * n_, n_square_) * (1 + s * n_) % n_square_;
    return Z;
}

LPuzzle LHTLP::GeneratePuzzle(const NTL::ZZ &s, const NTL::ZZ &r)
{
    assert(s < n_ && s >= 0);
    assert(r < n_ / 2 && r >= 0);
    LPuzzle Z(n_);
    Z.u = PowerMod(g_, r, n_);
    Z.v = PowerMod(h_, r * n_, n_square_) * (1 + s * n_) % n_square_;
    return Z;
}

// return solution s, will return -1 if the puzzle is invalid
NTL::ZZ LHTLP::QuickSolvePuzzle(const LPuzzle &Z)
{
    IS_CHEET_MODE
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &v = Z.v;
    NTL::ZZ w = NTL::PowerMod(u, trapdoor_, n_);
    NTL::ZZ temp = (v * InvMod(PowerMod(w, n_, n_square_), n_square_) % n_square_ - 1);
    return temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);
}

std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> LHTLP::QuickSolvePuzzleWithProof(const LPuzzle &Z)
{
    IS_CHEET_MODE
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &v = Z.v;
    NTL::ZZ w = NTL::PowerMod(u, trapdoor_, n_);
    NTL::ZZ temp = (v * InvMod(PowerMod(w, n_, n_square_), n_square_) % n_square_ - 1);
    NTL::ZZ s = temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);

    NTL::ZZ l = HashToPrime(u, w);

    NTL::ZZ r = PowerMod(NTL::ZZ(2), T_, l);
    NTL::ZZ q = (PowerMod(NTL::ZZ(2), T_, lambda_) - r) * InvMod(l, lambda_) % lambda_;
    NTL::ZZ pi = PowerMod(u, q, n_);
    return {s, w, pi};
}

NTL::ZZ LHTLP::SolvePuzzle(const LPuzzle &Z)
{
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &v = Z.v;
    NTL::ZZ w = u;
    for (NTL::ZZ i(0); i < T_; i++)
    {
        w = NTL::SqrMod(w, n_);
    }
    NTL::ZZ temp = (v * InvMod(PowerMod(w, n_, n_square_), n_square_) % n_square_ - 1);
    NTL::ZZ s = temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);
    return s;
}

std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> LHTLP::SolvePuzzleWithProof(const long k, const long gamma, const LPuzzle &Z)
{
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &v = Z.v;
    NTL::ZZ w = u;
    std::vector<NTL::ZZ> C(T_ / (k * gamma) + 1);
    for (long i = 0; i < T_; i++)
    {
        if (i % (k * gamma) == 0)
        {
            C[i / (k * gamma)] = w;
        }
        w = NTL::SqrMod(w, n_);
    }
    NTL::ZZ temp = v * PowerMod(w, -n_, n_square_) % n_square_ - 1;
    NTL::ZZ s = temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);

    NTL::ZZ l = HashToPrime(u, w);

    clock_t start_time = clock();
    NTL::ZZ pi = GenerateProof(k, gamma, C, l);
    clock_t end_time = clock(); 
    std::cout << "Time of generating a MCorSol/MIvalid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << std::endl;

    return {s, w, pi};
}

// invalid proof = 0, valid solution = 1, invalid solution = -1
int LHTLP::VerifyProofOfSol(const LPuzzle Z, const std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ> &solution_and_proof)
{
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &v = Z.v;
    const NTL::ZZ &s = std::get<0>(solution_and_proof);
    const NTL::ZZ &w = std::get<1>(solution_and_proof);
    const NTL::ZZ &pi = std::get<2>(solution_and_proof);

    NTL::ZZ l = HashToPrime(u, w);
    NTL::ZZ r = PowerMod(NTL::ZZ(2), T_, l);
    if (PowerMod(pi, l, n_) * PowerMod(u, r, n_) % n_ != w)
    {
        return 0;
    }
    if (s >= 0 && s < n_ && PowerMod(w, n_, n_square_) * (1 + s * n_) % n_square_ == v)
    {
        return 1;
    }

    if (s == -1 && (v * PowerMod(w, -n_, n_square_) % n_square_ - 1) % n_ != 0)
    {
        return -1;
    }

    return 0;
}
