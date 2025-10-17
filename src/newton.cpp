#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include "newton.h"

// Constructor with separate diagonals
NewtonSolver::NewtonSolver(
        const TridiagonalMatrix& sigma0Inv_val,
        const std::vector<double>& ks_val,
        const std::vector<double>& cs_val,
        const std::vector<double>& ms_val,
        double tol_val,
        int maxIters_val,
        double armijoC_val,
        double armijoRho_val,
        const std::vector<double>& qInit
    ) : 
    n(ks_val.size()),
    sigma0Inv(sigma0Inv_val), 
    ks(ks_val), 
    cs(cs_val),
    ms(ms_val),
    tol(tol_val),
    maxIters(maxIters_val),
    armijoC(armijoC_val),
    armijoRho(armijoRho_val),
    qCurr(qInit)
{
    // Initialize other variables to garbage values; they get updated
    eta = 0.0;
    gradCurr = std::vector<double>(n, tol + 1.0); // tol + 1 so that algortihm doesn't think it's already converged
    dCurr = std::vector<double>(n, 0.0);
    
    // To update all other private variables, take a step
    step();
}

double NewtonSolver::computeObjective(const std::vector<double> q) {

    // Remember this is the NEGATIVE log posterior 
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

    return result;
}

void NewtonSolver::computeGradient() {

    // 3. Add 1 + derivative(mj + log(kj + exp(-qj)))
    for(int i = 0; i < n; i++){
        gradCurr[i] = -((cs[i] * ms[i]) / (cs[i] + ks[i]*exp(qCurr[i])));

        if(i < n-1) {
            gradCurr[i]++;
            gradCurr[i] -= exp(qCurr[i+1] - qCurr[i]);
        }

        if(i > 0) {
            gradCurr[i] += exp(qCurr[i] - qCurr[i-1]);
        } 

        if(i == 0) {
            gradCurr[i]++;
        }
    }
}

TridiagonalMatrix NewtonSolver::computeHessian() {
    // 1. Start with sigma0Inv
    TridiagonalMatrix result = sigma0Inv;

    // sigma0Inv gets completely replaced

    // 2. Adjust diagonal
    for(int i = 0; i < n; i++){
        double denom = cs[i] + ks[i] * exp(qCurr[i]);
        double elem = ((cs[i] * ks[i] * ms[i] * exp(qCurr[i])) / (denom * denom));

        if(i < n-1) {
            elem += exp(qCurr[i+1] - qCurr[i]);
        }
        if(i > 0){
            elem += exp(qCurr[i] - qCurr[i-1]);
            // Set off diagonal too
            result.set(i, i-1, -exp(qCurr[i] - qCurr[i-1]));
            result.set(i-1, i, -exp(qCurr[i] - qCurr[i-1]));
        }
        result.set(i, i, elem);
    }

    return result;
}

void NewtonSolver::computeDescentDirection() {
    // Descent direction is -Hessian^-1 * gradient
    // Set y = -Hessian^-1 * gradient, then Hessian * y = -gradient and solve for y

    dCurr = computeHessian().solve(gradCurr);
    for(int i = 0; i < n; i++){
        dCurr[i] *= -1.0;
    }
}

// Compute ax + by, a, b scalars, x, y vectors
std::vector<double> vecSum(double a, const std::vector<double>& x, double b, const std::vector<double>& y) {
    std::vector<double> result = x;
    for(int i = 0; i < static_cast<int>(x.size()); i++) {
        result[i] *= a;
        result[i] += b * y[i];
    }
    return result;
}

double vecDot(const std::vector<double>& x, const std::vector<double>& y) {
    double result = 0.0;
    for(int i = 0; i < static_cast<int>(x.size()); i++) {
        result += x[i] * y[i];
    }
    return result;
}

void NewtonSolver::computeLearningRate() {

    // Max step size
    double alpha = 10.0; // Test if this is big enough?
    int nIters = 0;
    while(
        computeObjective(
            vecSum(1.0, qCurr, alpha, dCurr)
        )
        >
        computeObjective(qCurr) + armijoC * alpha * vecDot(gradCurr, dCurr)
    ) {
        alpha *= armijoRho;
        nIters++;
    }

    if(nIters == 0){
        throw std::runtime_error("Starting alpha was not big enough");
    }

    if(nIters == maxIters) {
        std::cout << "Warning: Armijo line search did not converge after " << maxIters << " iterations. Step size may not be optimal." << std::endl;
    }

    eta = alpha;

}

void NewtonSolver::updateQ() {
    qCurr = vecSum(1.0, qCurr, eta, dCurr);
}

void NewtonSolver::step() {
    // Initialie gradient
    computeGradient();

    // Initialize descent direction
    // Remember, descent direction always needs to get updated after gradient!
    computeDescentDirection();
    
    // Update eta using Armijo rule
    computeLearningRate();

    // Update qCurr
    updateQ();
}

double infinityNorm(const std::vector<double>& x){
    double result = 0.0;
    for(int i = 0; i < static_cast<int>(x.size()); i++) {
        result = std::max(result, std::abs(x[i]));
    }
    return result;
}

void NewtonSolver::minimize() {
    // Compute infinity norm of gradient
    while(infinityNorm(gradCurr) > tol && nSteps < maxIters) {
        step();
        nSteps++;
    }

    if(nSteps == maxIters) {
        std::cout << "Warning: Newton's method did not converge after " << maxIters << " iterations." << std::endl;
        std::cout << "Infinity norm of gradient: " << infinityNorm(gradCurr) << std::endl;
        std::cout << "Tolerance: " << tol << std::endl;
    }

    
}

std::vector<double> NewtonSolver::argmin() {
    return qCurr;
}

void NewtonSolver::print() {
    std::cout << "Number of iterations: " << nSteps << std::endl;
    std::cout << "Argmin:" << nSteps << std::endl;
    for(int i = 0; i < n; i++) {
        std::cout << qCurr[i] << std::endl;
    }
}