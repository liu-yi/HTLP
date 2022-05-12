#include "LHTLP.hpp"
#include "MHTLP.hpp"
#include "Puzzle.hpp"

#include <NTL/ZZ.h>
#include <vector>
#include <assert.h>
#include <openssl/sha.h>

using namespace std;
using namespace NTL;

clock_t start_time, end_time;

void TestLHTLP(const long k, const long gamma, const long T, const long MODULUS_LEN = 2048, const long kappa = 128)
{
    cout << "Init LHTLP scheme. " << endl;
    LHTLP scheme(MODULUS_LEN, T, kappa, true); // Init LHTLP scheme.

    ZZ s = scheme.GenerateRandomElement(); // Generate a random solution from Z_n
    ZZ s_prime = scheme.GenerateRandomElement();
    ZZ r = scheme.GenerateRandomExponent(); // Generate randomness for puzzle s
    LPuzzle Z = scheme.GeneratePuzzle(s);   // Generate a LHTLP for solution s
    LPuzzle Z_prime = scheme.GeneratePuzzle(s_prime);
    LPuzzle Z_LValid = scheme.GeneratePuzzle(s, r);

    LPuzzle Z_add = Z + Z_prime;

    LPuzzle invalid_Z(Z.n());
    invalid_Z.u = scheme.GenerateJacobiOne();
    invalid_Z.v = RandomBnd(scheme.n_square());

    auto sol_and_proof = scheme.SolvePuzzleWithProof(k, gamma, Z);
    auto invalid_Z_proof = scheme.SolvePuzzleWithProof(k, gamma, invalid_Z);
    auto sol_and_invalid_proof = sol_and_proof;
    get<2>(sol_and_invalid_proof) = scheme.GenerateRandomElement();

    auto s_add = scheme.SolvePuzzle(Z_add);
    assert((s + s_prime) % scheme.n() == s_add);

    auto LValid_proof = scheme.GenerateLValidProof(Z_LValid, s, r);

    // vector<double> time(3, 0);
    // for (int i = 0; i < 500; i++)
    // {
    start_time = clock();
    assert(scheme.VerifyProofOfSol(Z, sol_and_proof) == 1);
    end_time = clock();
    // time[0] += (double)(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Time of VerifyProofOfSol for 1 \t\t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

    start_time = clock();
    assert(scheme.VerifyProofOfSol(invalid_Z, invalid_Z_proof) == -1);
    end_time = clock();
    // time[1] += (double)(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Time of VerifyProofOfSol for -1 \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

    start_time = clock();
    assert(scheme.VerifyProofOfSol(Z, sol_and_invalid_proof) == 0);
    end_time = clock();
    // time[2] += (double)(end_time - start_time) / CLOCKS_PER_SEC;
    cout << "Time of VerifyProofOfSol for 0 \t\t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;
    // }
    // for (int i = 0; i < 3; i++)
    // {
    //     cout << i << " " << time[i] << endl;
    // }

    start_time = clock();
    assert(scheme.VerifyLValidProof(Z_LValid, LValid_proof) == true);
    end_time = clock();
    cout << "Time of VerifyLValidProof \t\t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;
}

void TestLCorSol(const long k, const long gamma, const long T, const long MODULUS_LEN = 2048, const long kappa = 128)
{
    cout << "Test LCorSol. " << endl;
    LHTLP scheme(MODULUS_LEN, T, kappa, true); // Init LHTLP scheme.

    ZZ s = scheme.GenerateRandomElement(); // Generate a random solution from Z_n
    LPuzzle Z = scheme.GeneratePuzzle(s);  // Generate a LHTLP for solution s

    start_time = clock();
    auto sol_and_proof = scheme.SolvePuzzleWithProof(k, gamma, Z);
    end_time = clock();
    cout << "Time of solving a puzzle and generating a LCorSol proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

    start_time = clock();
    assert(scheme.VerifyProofOfSol(Z, sol_and_proof) == 1);
    end_time = clock();
    cout << "Time of verifying a LCorSol proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;
}

void TestLInvalid(const long k, const long gamma, const long T, const long MODULUS_LEN = 2048, const long kappa = 128)
{
    cout << "Test LInvalid. " << endl;
    LHTLP scheme(MODULUS_LEN, T, kappa, true); // Init LHTLP scheme.

    ZZ s = scheme.GenerateRandomElement(); // Generate a random solution from Z_n
    LPuzzle Z = scheme.GeneratePuzzle(s);  // Generate a LHTLP for solution s

    LPuzzle invalid_Z(Z.n());
    invalid_Z.u = scheme.GenerateJacobiOne();
    invalid_Z.v = RandomBnd(scheme.n_square());

    start_time = clock();
    auto invalid_Z_proof = scheme.SolvePuzzleWithProof(k, gamma, invalid_Z);
    end_time = clock();
    cout << "Time of solving a puzzle and generating a LInvalid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

    start_time = clock();
    assert(scheme.VerifyProofOfSol(invalid_Z, invalid_Z_proof) == -1);
    end_time = clock();
    cout << "Time of verifying a LInvalid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;
}

void TestLValid(const long k, const long gamma, const long T, const long MODULUS_LEN = 2048, const long kappa = 128)
{
    cout << "Test LValid. " << endl;
    LHTLP scheme(MODULUS_LEN, T, kappa, true); // Init LHTLP scheme.

    ZZ s = scheme.GenerateRandomElement();  // Generate a random solution from Z_n
    ZZ r = scheme.GenerateRandomExponent(); // Generate randomness for puzzle s
    LPuzzle Z_LValid = scheme.GeneratePuzzle(s, r);

    start_time = clock();
    auto LValid_proof = scheme.GenerateLValidProof(Z_LValid, s, r);
    end_time = clock();
    cout << "Time of generating a LValid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

    start_time = clock();
    assert(scheme.VerifyLValidProof(Z_LValid, LValid_proof) == true);
    end_time = clock();
    cout << "Time of verifying a LValid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;
}

void TestMHTLP(const long k, const long gamma, const long T, const long MODULUS_LEN = 2048, const long kappa = 128)
{
    cout << "Init MHTLP scheme. " << endl;
    MHTLP scheme(MODULUS_LEN, T, kappa, true); // Init LHTLP scheme.

    ZZ s = scheme.GenerateRandomElement(); // Generate a random solution from Z_n
    ZZ s_prime = scheme.GenerateRandomElement();
    ZZ r = scheme.GenerateRandomExponent();       // Custom randomness
    ZZ r_prime = scheme.GenerateRandomExponent(); // Custom randomness
    MPuzzle Z = scheme.GeneratePuzzle(s);         // Generate a LHTLP for solution s
    MPuzzle Z_prime = scheme.GeneratePuzzle(s_prime);
    MPuzzle Z_MValid = scheme.GeneratePuzzle(s, r, r_prime);

    MPuzzle Z_mul = Z * Z_prime;

    MPuzzle invalid_Z(Z.n());
    invalid_Z.u = scheme.GenerateJacobiOne();
    invalid_Z.u_prime = scheme.GenerateJacobiOne();
    invalid_Z.v = scheme.GenerateJacobiOne();
    invalid_Z.theta = RandomBnd(scheme.n_square());

    // cout << "hello" << endl;
    auto sol_and_proof = scheme.SolvePuzzleWithProof(k, gamma, Z);
    // auto sol_and_proof = scheme.QuickSolvePuzzleWithProof(Z);
    auto invalid_Z_proof = scheme.SolvePuzzleWithProof(k, gamma, invalid_Z);
    auto sol_and_invalid_proof = sol_and_proof;
    get<2>(sol_and_invalid_proof) = scheme.GenerateRandomElement();

    auto s_mul = scheme.SolvePuzzle(Z_mul);
    assert(s * s_prime % scheme.n() == s_mul);

    auto MValid_proof = scheme.GenerateMValidProof(Z_MValid, s, r, r_prime);

    assert(s == std::get<0>(sol_and_proof));
    assert(-1 == std::get<0>(invalid_Z_proof));

    assert(scheme.VerifyProofOfSol(Z, sol_and_proof) == 1);
    assert(scheme.VerifyProofOfSol(invalid_Z, invalid_Z_proof) == -1);
    assert(scheme.VerifyProofOfSol(Z, sol_and_invalid_proof) == 0);

    assert(scheme.VerifyMValidProof(Z_MValid, MValid_proof) == true);
}

void TestMCorSol(const long k, const long gamma, const long T, const long MODULUS_LEN = 2048, const long kappa = 128)
{
    cout << "Init MHTLP scheme. " << endl;
    MHTLP scheme(MODULUS_LEN, T, kappa, true); // Init LHTLP scheme.

    ZZ s = scheme.GenerateRandomElement(); // Generate a random solution from Z_n
    MPuzzle Z = scheme.GeneratePuzzle(s);  // Generate a LHTLP for solution s

    start_time = clock();
    auto sol_and_proof = scheme.SolvePuzzleWithProof(k, gamma, Z);
    end_time = clock();
    cout << "Time of solving a puzzle and generating a MCorSol proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

    start_time = clock();
    assert(scheme.VerifyProofOfSol(Z, sol_and_proof) == 1);
    end_time = clock();
    cout << "Time of verifying a MCorSol proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;
}

void TestMInvalid(const long k, const long gamma, const long T, const long MODULUS_LEN = 2048, const long kappa = 128)
{
    cout << "Test MInvalid. " << endl;
    MHTLP scheme(MODULUS_LEN, T, kappa, true); // Init LHTLP scheme.

    ZZ s = scheme.GenerateRandomElement(); // Generate a random solution from Z_n
    MPuzzle Z = scheme.GeneratePuzzle(s);  // Generate a LHTLP for solution s

    MPuzzle invalid_Z(Z.n());
    invalid_Z.u = scheme.GenerateJacobiOne();
    invalid_Z.u_prime = scheme.GenerateJacobiOne();
    invalid_Z.v = scheme.GenerateJacobiOne();
    invalid_Z.theta = RandomBnd(scheme.n_square());

    start_time = clock();
    auto invalid_Z_proof = scheme.SolvePuzzleWithProof(k, gamma, invalid_Z);
    end_time = clock();
    cout << "Time of solving a puzzle and generating a MInvalid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

    start_time = clock();
    assert(scheme.VerifyProofOfSol(invalid_Z, invalid_Z_proof) == -1);
    end_time = clock();
    cout << "Time of verifying a MInvalid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;
}

void TestMValid(const long k, const long gamma, const long T, const long MODULUS_LEN = 2048, const long kappa = 128)
{
    cout << "Test MValid. " << endl;
    MHTLP scheme(MODULUS_LEN, T, kappa, true); // Init LHTLP scheme.

    ZZ s = scheme.GenerateRandomElement();        // Generate a random solution from Z_n
    ZZ r = scheme.GenerateRandomExponent();       // Custom randomness
    ZZ r_prime = scheme.GenerateRandomExponent(); // Custom randomness
    MPuzzle Z_MValid = scheme.GeneratePuzzle(s, r, r_prime);

    start_time = clock();
    auto MValid_proof = scheme.GenerateMValidProof(Z_MValid, s, r, r_prime);
    end_time = clock();
    cout << "Time of solving a puzzle and generating a MValid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;

    start_time = clock();
    assert(scheme.VerifyMValidProof(Z_MValid, MValid_proof) == true);
    end_time = clock();
    cout << "Time of verifying a MValid proof \t" << (double)(end_time - start_time) / CLOCKS_PER_SEC << endl;
}

int main(void)
{
    // Parameter setting for VDF
    const long T = 10000000;
    const long k = max(1, long(log2(T)) / 3);
    const long gamma = max(1, long(sqrt(T)));

    // TestLCorSol(k, gamma, T);

    // TestLInvalid(k, gamma, T);

    // TestLValid(k, gamma, T);

    // TestMCorSol(k, gamma, T);

    // TestMInvalid(k, gamma, T);

    // TestMValid(k, gamma, T);

    // TestLHTLP(k, gamma, T);

    // TestMHTLP(k, gamma, T);

    return 0;
}