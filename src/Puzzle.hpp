#ifndef PUZZLE_HPP_
#define PUZZLE_HPP_


class APuzzle
{
private:
    NTL::ZZ n_;

public:
    NTL::ZZ u;
    NTL::ZZ v;
    APuzzle(const NTL::ZZ &u, const NTL::ZZ &v, const NTL::ZZ n) : n_(n)
    {
        this->u = u;
        this->v = v;
    }
    APuzzle(const NTL::ZZ n) : n_(n)
    {
    }

    NTL::ZZ n()
    {
        return n_;
    }

    APuzzle operator+(const APuzzle &a)
    {
        assert(n_ == a.n_);
        APuzzle Z(u * a.u % n_, v * a.v % (n_ * n_), n_);
        return Z;
    }

    friend APuzzle operator+(const NTL::ZZ &a, const APuzzle &b)
    {
        NTL::ZZ c = a % b.n_;
        APuzzle Z(b.u, (1 + c * b.n_) * b.v % (b.n_ * b.n_), b.n_);
        return Z;
    }

    friend APuzzle operator+(long a, const APuzzle &b)
    {
        NTL::ZZ c = NTL::ZZ(a) % b.n_;
        APuzzle Z(b.u, (1 + c * b.n_) * b.v % (b.n_ * b.n_), b.n_);
        return Z;
    }

    friend APuzzle operator*(long a, const APuzzle &b)
    {
        APuzzle Z(PowerMod(b.u, a, b.n_), PowerMod(b.v, a, b.n_) % (b.n_ * b.n_), b.n_);
        return Z;
    }
};

class MPuzzle
{
private:
    NTL::ZZ n_;

public:
    NTL::ZZ u;
    NTL::ZZ u_prime;
    NTL::ZZ v;
    NTL::ZZ theta;

    MPuzzle(const NTL::ZZ &u, const NTL::ZZ &u_prime, const NTL::ZZ &v, const NTL::ZZ &theta, const NTL::ZZ n) : n_(n)
    {
        this->u = u;
        this->u_prime = u_prime;
        this->v = v;
        this->theta = theta;
    }
    MPuzzle(const NTL::ZZ n) : n_(n)
    {
    }

    MPuzzle operator*(const MPuzzle &a)
    {
        assert(n_ == a.n_);
        MPuzzle Z(u * a.u % n_, u_prime * a.u_prime % n_, v * a.v % n_, theta * a.theta % (n_ * n_), n_);
        return Z;
    }

    MPuzzle Power(const long e){
        MPuzzle Z(PowerMod(u, e, n_), PowerMod(u_prime, e, n_), PowerMod(v, e, n_), PowerMod(theta, e, n_), n_); 
        return Z; 
    }

    NTL::ZZ n()
    {
        return n_;
    }
};

#endif