#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include "matrix.h"

class NewtonSolver {
private:

    const int n; // Noting this will usually be n-1 in our context; it's simply the dimension of the space here

    const TridiagonalMatrix sigma0Inv;
    const std::vector<double> ks;
    const std::vector<double> cs;
    const std::vector<double> ms;

    // Tolerance
    const double tol;
    const double maxIters;
    // Armijo params
    const double armijoC;
    const double armijoRho;
    
    // Current value of iterate
    std::vector<double> qCurr;

    // Learning rate
    double eta;

    // Current value of gradient
    std::vector<double> gradCurr;

    // Current value of descent direction
    std::vector<double> dCurr;

    // Number of iterations
    int nSteps = 0;
    

public:

    // Constructor with separate diagonals
    NewtonSolver(
        const TridiagonalMatrix& sigma0Inv_val,
        const std::vector<double>& ks_val,
        const std::vector<double>& cs_val,
        const std::vector<double>& ms_val,
        double tol_val,
        int maxIters_val,
        double armijoC_val,
        double armijoRho_val,
        const std::vector<double>& qInit
    );

    double computeObjective(const std::vector<double> q);
    void computeGradient(); // Evaluated at qCurr
    TridiagonalMatrix computeHessian(); // Evaluated at qCurr
    void computeDescentDirection(); // Evaluated at qCurr
    void computeLearningRate(); // Updates eta
    void updateQ(); // Updates qCurr

    void step(); // One full iteration

    void minimize();

    std::vector<double> argmin();

    void print();
};
