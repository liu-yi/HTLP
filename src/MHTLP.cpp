#include "MHTLP.hpp"

#define IS_CHEET_MODE assert(cheeting_mode_);

MHTLP::MHTLP(const long modulus_len, const long T, const long kappa) : HTLP(modulus_len, T, kappa)
{
    while (true)
    {
        chi_ = GenerateRandomElement();
        if (Jacobi(chi_, n_) == -1)
        {
            break;
        }
    }
}
MHTLP::MHTLP(const NTL::ZZ &n, const NTL::ZZ &g, const NTL::ZZ &h, const NTL::ZZ &chi, const long T, const long kappa) : HTLP(n, g, h, T, kappa)
{
    chi_ = chi;
}
MHTLP::MHTLP(const long modulus_len, const long T, const long kappa, bool cheeting_mode) : HTLP(modulus_len, T, kappa, cheeting_mode)
{
    while (true)
    {
        chi_ = GenerateRandomElement();
        if (Jacobi(chi_, n_) == -1)
        {
            break;
        }
    }
}

std::tuple<std::vector<NTL::ZZ>, std::vector<NTL::ZZ>, std::vector<NTL::ZZ>, std::vector<NTL::ZZ>> MHTLP::GenerateMValidProof(const MPuzzle &Z, const NTL::ZZ &s, const NTL::ZZ &r, const NTL::ZZ &r_prime)
{
    int sigma;
    if (NTL::Jacobi(s, n_) == -1)
    {
        sigma = 1;
    }
    else
    {
        sigma = 0;
    }
    std::vector<NTL::ZZ> theta(2);
    theta[0] = Z.theta;
    theta[1] = Z.theta * NTL::InvMod(1 + n_, n_square_) % n_square_;
    const NTL::ZZ &u_prime = Z.u_prime;
    std::vector<NTL::ZZ> a(2);
    std::vector<NTL::ZZ> b(2);
    std::vector<NTL::ZZ> e(2);
    std::vector<NTL::ZZ> alpha(2);
    std::vector<NTL::ZZ> beta(2);

    e[1 - sigma] = NTL::RandomLen_ZZ(kappa_);
    alpha[1 - sigma] = NTL::RandomLen_ZZ(modulus_len_ - 1 + 2 * kappa_);
    a[1 - sigma] = PowerMod(g_, alpha[1 - sigma], n_) * PowerMod(u_prime, -e[1 - sigma], n_) % n_;
    b[1 - sigma] = PowerMod(h_, alpha[1 - sigma] * n_, n_square_) * PowerMod(theta[1 - sigma], -e[1 - sigma], n_square_) % n_square_;

    NTL::ZZ x = NTL::RandomLen_ZZ(modulus_len_ - 1 + 2 * kappa_);
    a[sigma] = PowerMod(g_, x, n_);
    b[sigma] = PowerMod(h_, x * n_, n_square_);

    std::stringstream ss;
    ss << n_ << g_ << h_ << Z.u << Z.u_prime << Z.v << Z.theta << a[0] << a[1] << b[0] << b[1];
    NTL::ZZ e_challenge = HashToElement(ss.str());

    e[sigma] = e_challenge ^ e[1 - sigma];

    alpha[sigma] = r_prime * e[sigma] + x;

    return {a, b, e, alpha};
}

bool MHTLP::VerifyMValidProof(const MPuzzle &Z, const std::tuple<std::vector<NTL::ZZ>, std::vector<NTL::ZZ>, std::vector<NTL::ZZ>, std::vector<NTL::ZZ>> &proof)
{
    const std::vector<NTL::ZZ> &a = std::get<0>(proof);
    const std::vector<NTL::ZZ> &b = std::get<1>(proof);
    const std::vector<NTL::ZZ> &e = std::get<2>(proof);
    const std::vector<NTL::ZZ> &alpha = std::get<3>(proof);
    std::vector<NTL::ZZ> theta(2);
    theta[0] = Z.theta;
    theta[1] = Z.theta * NTL::InvMod(1 + n_, n_square_) % n_square_;
    std::stringstream ss;
    ss << n_ << g_ << h_ << Z.u << Z.u_prime << Z.v << Z.theta << a[0] << a[1] << b[0] << b[1];
    NTL::ZZ e_challenge = HashToElement(ss.str());
    if (e_challenge != (e[0] ^ e[1]))
    {

        return false;
    }
    for (int i = 0; i < 2; i++)
    {
        if (PowerMod(Z.u_prime, e[i], n_) * a[i] % n_ != PowerMod(g_, alpha[i], n_))
        {
            return false;
        }
        if (PowerMod(theta[i], e[i], n_square_) * b[i] % n_square_ != PowerMod(h_, alpha[i] * n_, n_square_))
        {
            return false;
        }
    }
    return true;
}

MPuzzle MHTLP::GeneratePuzzle(const NTL::ZZ &s)
{
    assert(s < n_ && s >= 0);
    MPuzzle Z(n_);
    NTL::ZZ r = GenerateRandomExponent();
    NTL::ZZ r_prime = GenerateRandomExponent();
    int sigma;
    if (NTL::Jacobi(s, n_) == -1)
    {
        sigma = 1;
    }
    else
    {
        sigma = 0;
    }
    Z.u = PowerMod(g_, r, n_);
    Z.u_prime = PowerMod(g_, r_prime, n_);
    Z.v = PowerMod(h_, r, n_) * PowerMod(chi_, sigma, n_) * s % n_;
    Z.theta = PowerMod(h_, r_prime * n_, n_square_) * (1 + sigma * n_) % n_square_;
    return Z;
}
MPuzzle MHTLP::GeneratePuzzle(const NTL::ZZ &s, const NTL::ZZ &r, const NTL::ZZ &r_prime)
{
    assert(s < n_ && s >= 0);
    assert(r < n_ / 2 && r >= 0);
    assert(r_prime < n_ / 2 && r_prime >= 0);
    MPuzzle Z(n_);
    int sigma;
    if (NTL::Jacobi(s, n_) == -1)
    {
        sigma = 1;
    }
    else
    {
        sigma = 0;
    }
    Z.u = PowerMod(g_, r, n_);
    Z.u_prime = PowerMod(g_, r_prime, n_);
    Z.v = PowerMod(h_, r, n_) * PowerMod(chi_, sigma, n_) * s % n_;
    Z.theta = PowerMod(h_, r_prime * n_, n_square_) * (1 + sigma * n_) % n_square_;
    return Z;
}

NTL::ZZ MHTLP::SolvePuzzle(const MPuzzle &Z)
{
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &u_prime = Z.u_prime;
    const NTL::ZZ &v = Z.v;
    const NTL::ZZ &theta = Z.theta;
    NTL::ZZ w = u;
    NTL::ZZ w_prime = u_prime;
    for (NTL::ZZ i(0); i < T_; i++)
    {
        w = NTL::SqrMod(w, n_);
    }
    for (NTL::ZZ i(0); i < T_; i++)
    {
        w_prime = NTL::SqrMod(w_prime, n_);
    }
    NTL::ZZ temp = theta * PowerMod(w_prime, -n_, n_square_) % n_square_ - 1;
    NTL::ZZ sigma = temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);

    NTL::ZZ s;
    if (sigma == -1)
    {
        s = -1;
    }
    else
    {
        s = v * PowerMod(chi_, -sigma, n_) * NTL::InvMod(w, n_) % n_;
    }
    return s;
}

// return solution s, will return -1 if the puzzle is invalid
NTL::ZZ MHTLP::QuickSolvePuzzle(const MPuzzle &Z)
{
    IS_CHEET_MODE
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &u_prime = Z.u_prime;
    const NTL::ZZ &v = Z.v;
    const NTL::ZZ &theta = Z.theta;

    NTL::ZZ w = NTL::PowerMod(u, trapdoor_, n_);
    NTL::ZZ w_prime = NTL::PowerMod(u_prime, trapdoor_, n_);
    NTL::ZZ temp = theta * PowerMod(w_prime, -n_, n_square_) % n_square_ - 1;
    NTL::ZZ sigma = temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);

    NTL::ZZ s;
    if (sigma == -1)
    {
        s = -1;
    }
    else
    {
        s = v * PowerMod(chi_, -sigma, n_) * NTL::InvMod(w, n_) % n_;
    }
    return s;
}

std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> MHTLP::QuickSolvePuzzleWithProof(const MPuzzle &Z)
{
    IS_CHEET_MODE
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &u_prime = Z.u_prime;
    const NTL::ZZ &v = Z.v;
    const NTL::ZZ &theta = Z.theta;

    NTL::ZZ w = NTL::PowerMod(u, trapdoor_, n_);
    NTL::ZZ w_prime = NTL::PowerMod(u_prime, trapdoor_, n_);
    NTL::ZZ temp = theta * PowerMod(w_prime, -n_, n_square_) % n_square_ - 1;
    NTL::ZZ sigma = temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);

    NTL::ZZ s;
    if (sigma == -1)
    {
        s = -1;
    }
    else
    {
        s = v * PowerMod(chi_, -sigma, n_) * NTL::InvMod(w, n_) % n_;
    }

    NTL::ZZ l = HashToPrime(u, w);
    NTL::ZZ l_prime = HashToPrime(u_prime, w_prime);

    NTL::ZZ r = PowerMod(NTL::ZZ(2), T_, l);
    NTL::ZZ r_prime = PowerMod(NTL::ZZ(2), T_, l_prime);
    NTL::ZZ q = (PowerMod(NTL::ZZ(2), T_, lambda_) - r) * InvMod(l, lambda_) % lambda_;
    NTL::ZZ q_prime = (PowerMod(NTL::ZZ(2), T_, lambda_) - r_prime) * InvMod(l_prime, lambda_) % lambda_;
    NTL::ZZ pi = PowerMod(u, q, n_);
    NTL::ZZ pi_prime = PowerMod(u_prime, q_prime, n_);
    return {s, w, w_prime, pi, pi_prime};
}

std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> MHTLP::SolvePuzzleWithProof(const long k, const long gamma, const MPuzzle &Z)
{
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &u_prime = Z.u_prime;
    const NTL::ZZ &v = Z.v;
    const NTL::ZZ &theta = Z.theta;
    NTL::ZZ w = u;
    NTL::ZZ w_prime = u_prime;
    std::vector<NTL::ZZ> C(T_ / (k * gamma) + 1);
    std::vector<NTL::ZZ> C_prime(T_ / (k * gamma) + 1);
    for (long i = 0; i < T_; i++)
    {
        if (i % (k * gamma) == 0)
        {
            C[i / (k * gamma)] = w;
        }
        w = NTL::SqrMod(w, n_);
    }
    clock_t start = clock();
    for (long i = 0; i < T_; i++)
    {
        if (i % (k * gamma) == 0)
        {
            C_prime[i / (k * gamma)] = w_prime;
        }
        w_prime = NTL::SqrMod(w_prime, n_);
    }
        clock_t end = clock(); 
    std::cout << "Time  \t" << (double)(end - start) / CLOCKS_PER_SEC << std::endl;


    NTL::ZZ temp = theta * PowerMod(w_prime, -n_, n_square_) % n_square_ - 1;
    NTL::ZZ sigma = temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);

    NTL::ZZ s;
    if (sigma == -1)
    {
        s = -1;
    }
    else
    {
        s = v * PowerMod(chi_, -sigma, n_) * NTL::InvMod(w, n_) % n_;
    }

    NTL::ZZ l = HashToPrime(u, w);
    NTL::ZZ l_prime = HashToPrime(u_prime, w_prime);

    clock_t start_time = clock();

    NTL::ZZ pi = GenerateProof(k, gamma, C, l);
    NTL::ZZ pi_prime = GenerateProof(k, gamma, C_prime, l_prime);

    clock_t end_time = clock(); 
    std::cout << "Time of generating a MCorSol/MIvalid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << std::endl;

    return {s, w, w_prime, pi, pi_prime};
}

// invalid proof = 0, valid solution = 1, invalid solution = -1
int MHTLP::VerifyProofOfSol(const MPuzzle Z, const std::tuple<NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ, NTL::ZZ> &solution_and_proof)
{
    const NTL::ZZ &u = Z.u;
    const NTL::ZZ &u_prime = Z.u_prime;
    const NTL::ZZ &v = Z.v;
    const NTL::ZZ &theta = Z.theta;
    const NTL::ZZ &s = std::get<0>(solution_and_proof);
    const NTL::ZZ &w = std::get<1>(solution_and_proof);
    const NTL::ZZ &w_prime = std::get<2>(solution_and_proof);
    const NTL::ZZ &pi = std::get<3>(solution_and_proof);
    const NTL::ZZ &pi_prime = std::get<4>(solution_and_proof);

    NTL::ZZ l = HashToPrime(u, w);
    NTL::ZZ l_prime = HashToPrime(u_prime, w_prime);
    NTL::ZZ r = PowerMod(NTL::ZZ(2), T_, l);
    NTL::ZZ r_prime = PowerMod(NTL::ZZ(2), T_, l_prime);
    if (PowerMod(pi, l, n_) * PowerMod(u, r, n_) % n_ != w || PowerMod(pi_prime, l_prime, n_) * PowerMod(u_prime, r_prime, n_) % n_ != w_prime)
    {
        return 0;
    }

    NTL::ZZ temp = theta * PowerMod(w_prime, -n_, n_square_) % n_square_ - 1;
    NTL::ZZ sigma = temp % n_ == NTL::ZZ(0) ? temp / n_ : NTL::ZZ(-1);
    if (sigma == -1 && s == -1)
    {
        return -1;
    }
    if (s >= 0 && s <= n_ && v == w * PowerMod(chi_, sigma, n_) * s % n_)
    {
        return 1;
    }
    return 0;
}


