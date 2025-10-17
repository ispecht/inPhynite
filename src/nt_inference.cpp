#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include "matrix.h"
#include "tree_inference.h"
#include "newton.h"
#include "nt_inference.h"


void updateNt(
    State& state,
    std::mt19937& rng
) {

    // NUMBER OF LEAVES
    int n = state.n;

    // m[i] := number of intervals on which the effective population size is nt[i]
    std::vector<double> m(n-1, 0.0);

    for(size_t i = 0; i < state.k.size(); i++) {
        int intervalIdx = n - state.k[i];
        m[intervalIdx] += 1.0;
    }

    
    std::vector<double> k(n-1, 0.0);
    std::vector<double> c(n-1, 0.0);

    for(int j = 0; j < n-1; j++){
        k[j] = (double) n - j;
        c[j] = (k[j] * (k[j] - 1.0)) / 2.0;
    }

    // sInv is just any matrix of the right size
    TridiagonalMatrix sInv(
        std::vector<double>(n-1, 0.0),
        std::vector<double>(n-2, 0.0)
    );

    // Initialize at current value
    std::vector<double> qInit = state.nt;
    for(size_t j = 0; j < qInit.size(); j++) {
        qInit[j] = std::log(qInit[j]);
    }

    // Newton's method
    NewtonSolver newton(
        sInv,
        k,
        c,
        m,
        0.0001,
        10000,
        0.0001,
        0.5,
        qInit
    );

    newton.minimize();


    // 1. Obtain approximate mean by way of Newton's method
    std::vector<double> b = newton.argmin();

    // Completely replaces sInv
    for(int i = 0; i < n-1; i++){
        double denom = c[i] + k[i] * exp(b[i]);
        double elem = ((c[i] * k[i] * m[i] * exp(b[i])) / (denom * denom));

        if(i < n-1) {
            elem += exp(b[i+1] - b[i]);
        }
        if(i > 0){
            elem += exp(b[i] - b[i-1]);
            // Set off diagonal too
            sInv.set(i, i-1, -exp(b[i] - b[i-1]));
            sInv.set(i-1, i, -exp(b[i] - b[i-1]));
        }
        sInv.set(i, i, elem);
    }


    // For sampling, need Cholesky decomp of sInv
    std::pair<std::vector<double>, std::vector<double>> L = sInv.choleskyDecomposition();

    // Have S = L^-T L^-1. Hence want to compute L^-T z, z[j] standard normal.
    // To do this set y = L^-T z, solve L^T y = z for y
    std::normal_distribution<double> normal_dist(0.0, 1.0); // Mean=0, StdDev=1
    
    std::vector<double> z(n - 1);
    
    // Generate n-1 standard normal random variables
    for (int j = 0; j < n - 1; ++j) {
        z[j] = normal_dist(rng);
    }

    // Get y = L^-T z, i.e. y solves L^T y = z
    std::vector<double> y = solveUpperBidiagonal(L, z);

    // Add the mean
    for(int j = 0; j < n - 1; j++){
        y[j] += b[j];
    }

    // Exponentiate
    for(int j = 0; j < n - 1; j++){
        y[j] = exp(y[j]);
    }

    state.nt = y;

}