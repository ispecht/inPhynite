#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cmath>

class TridiagonalMatrix {
private:
    std::vector<double> lower_diag;  // Lower diagonal (n-1 elements)
    std::vector<double> main_diag;   // Main diagonal (n elements)
    std::vector<double> upper_diag;  // Upper diagonal (n-1 elements)
    int n;  // Matrix dimension

public:

    // Constructor with separate diagonals
    TridiagonalMatrix(const std::vector<double>& lower,
                        const std::vector<double>& main,
                        const std::vector<double>& upper);

    // Constructor for symmetric tridiagonal (off-diagonal elements are same)
    TridiagonalMatrix(const std::vector<double>& main,
                        const std::vector<double>& off_diag);

    // Constructor for scalar tridiagonal (constant values)
    TridiagonalMatrix(int size, double lower_val, double main_val, double upper_val);

    // Copy constructor
    TridiagonalMatrix(const TridiagonalMatrix& other);

    // Getters
    int size() const;
    const std::vector<double>& getLowerDiagonal() const;
    const std::vector<double>& getMainDiagonal() const;
    const std::vector<double>& getUpperDiagonal() const;

    // Element access (i, j)
    double operator()(int i, int j) const;

    void set(int i, int j, double value);

    // Solve Ax = b using Thomas algorithm
    std::vector<double> solve(const std::vector<double>& b) const;

    // Matrix-vector multiplication: y = Ax
    std::vector<double> multiply(const std::vector<double>& x) const;

    // Check if matrix is symmetric
    bool isSymmetric() const;

    bool isPositiveDefinite() const;

    std::pair<std::vector<double>, std::vector<double>> choleskyDecomposition() const;


    // Solve using Cholesky decomposition: A = L * L^T, solve L * y = b, then L^T * x = y
    std::vector<double> solveCholesky(const std::vector<double>& b) const;

    

    // Print matrix (for debugging)
    void print(std::ostream& os = std::cout) const;
    
};

std::vector<double> solveLowerBidiagonal(
    const std::pair<std::vector<double>, std::vector<double>>& L,
    const std::vector<double>& b);

std::vector<double> solveUpperBidiagonal(
    const std::pair<std::vector<double>, std::vector<double>>& U,
    const std::vector<double>& b);

void demonstrateTridiagonalMatrix();

