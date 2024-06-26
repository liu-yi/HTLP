﻿# Homomorphic Time-Lock Puzzle Schemes 
![](https://img.shields.io/badge/version-v0.1-blue)

This is the repository for the paper 

Towards Practical Homomorphic Time-Lock Puzzles: Applicability and Verifiability, *ESORICS 2022* [[Link](https://link.springer.com/chapter/10.1007/978-3-031-17140-6_21)]

by Yi Liu, Qi Wang, Siu-Ming Yiu. 

## Introduction

This repository provides implementations for the **additively** homomorphic time-lock puzzle scheme with solution space $\mathbb{Z}_n$ in [[MT19]](https://eprint.iacr.org/2019/635.pdf) and our **multiplicatively** homomorphic time-lock puzzle scheme with solution space $\mathbb{Z}_n^*$. 

To avoid the redundant cost of the puzzle-solving process, we provide three *simple* and *fast* protocols for both the additively HTLP scheme with the solution space $\mathbb{Z}_n$ and our multiplicatively HTLP scheme, respectively, to verify the following three properties.

1. **Correctness.** A puzzle solver is able to convince other parties of the *correctness of the solution* that he solves from a puzzle.
2. **Invalidity.** Upon finding that a puzzle is invalid, one can convince other parties of the *invalidity of the puzzle*.
3. **Validity.** A puzzle generator can convince other parties of the *validity of the puzzle* he generated. 

## Dependencies

This project has dependencies of [NTL](https://github.com/libntl/ntl) and [OpenSSL](https://www.openssl.org/). 

## Build

1. Install [NTL](https://github.com/libntl/ntl). 
2. Install libssl-dev
3. Clone the repository: 
    ```
    git clone https://github.com/liu-yi/HTLP
    ```
4. Enter the directory
   ```
   cd HTLP
   ```
5. Assuming you have globally installed NTL and libssl-dev:
   ```
   make
   ```

6. You can test the execution by `./TestHTLP`, which is defined in `src/main.cpp`. 

