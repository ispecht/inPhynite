#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include "matrix.h"
#include "tree_inference.h"
#include "newton.h"
#include "nt_inference.h"

double getLogPosterior(
    const std::vector<double> q,
    const std::vector<double> ms,
    const std::vector<double> cs,
    const std::vector<double> ks
){
    int n = q.size();
    double result = 0.0;

    for(int i = 0; i < n; i++){
        result += q[i];
        result += ms[i] * log(ks[i] + cs[i] * exp(-q[i]));

        if(i >= 1) {
            result -= q[i];
            result += q[i-1];
            result += exp(q[i] - q[i-1]);
        }
    }



    return (-1.0 * result);
}

void updateNt(
    State& state,
    std::mt19937& rng,
    TridiagonalMatrix& sInv,
    bool forceAccept
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
    /*
    TridiagonalMatrix sInv(
        std::vector<double>(n-1, 0.0),
        std::vector<double>(n-2, 0.0)
    );
    */

    // Initialize at current value
    std::vector<double> qInit = state.nt;
    for(size_t j = 0; j < qInit.size(); j++) {
        qInit[j] = std::log(qInit[j]);
    }

    double log_posterior_old = getLogPosterior(qInit, m, c, k);

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

    // Current value of normal proposals z
    std::vector<double> z_curr = qInit;
    // Subtract the mean
    for(int j = 0; j < n - 1; j++){
        z_curr[j] -= b[j];
    }

    // have z_curr = L^T(qInit - b), solve 
    for(int j = 0; j < n-2; j++){
        z_curr[j] = (z_curr[j] * L.first[j]) + (z_curr[j+1] * L.second[j]);
    }
    z_curr[n-2] *= L.first[n-2];

    // Log probability of this proposal
    double log_p_new_to_old = 0.0;
    for(int j = 0; j < n-1; j++) {
        log_p_new_to_old -= (0.5 * z_curr[j] * z_curr[j]);
    }

    // Have S = L^-T L^-1. Hence want to compute L^-T z, z[j] standard normal.
    // To do this set y = L^-T z, solve L^T y = z for y
    std::normal_distribution<double> normal_dist(0.0, 1.0); // Mean=0, StdDev=1
    
    std::vector<double> z(n - 1);

    double log_p_old_to_new = 0.0;
    
    // Generate n-1 standard normal random variables
    for (int j = 0; j < n - 1; ++j) {
        z[j] = normal_dist(rng);
        log_p_old_to_new -= (0.5 * z[j] * z[j]);
    }

    // Get y = L^-T z, i.e. y solves L^T y = z
    std::vector<double> y = solveUpperBidiagonal(L, z);

    // Add the mean
    for(int j = 0; j < n - 1; j++){
        y[j] += b[j];
    }

    double log_posterior_new = getLogPosterior(y, m, c, k);

    // MH Accept/Reject

    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double u = uniform(rng);

    double log_p_accept = log_posterior_new - log_posterior_old + log_p_new_to_old - log_p_old_to_new;

    /*
    std::cout << "log_posterior_new: " << log_posterior_new << std::endl;
    std::cout << "log_posterior_old: " << log_posterior_old << std::endl;
    std::cout << "log_p_new_to_old: " << log_p_new_to_old << std::endl;
    std::cout << "log_p_old_to_new: " << log_p_old_to_new << std::endl;
    std::cout << "log_p_accept: " << log_p_accept << std::endl;
    */
    

    if((std::log(u) < log_p_accept) || forceAccept) {
        // Exponentiate
        for(int j = 0; j < n - 1; j++){
            y[j] = exp(y[j]);
        }
        state.nt = y; // Update state
    }
}