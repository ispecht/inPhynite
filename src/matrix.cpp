#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include "matrix.h"


// Constructor with separate diagonals
TridiagonalMatrix::TridiagonalMatrix(const std::vector<double>& lower,
                    const std::vector<double>& main,
                    const std::vector<double>& upper)
    : lower_diag(lower), main_diag(main), upper_diag(upper), n(main.size()) {
    
    if (static_cast<int>(lower.size()) != n - 1 || static_cast<int>(upper.size()) != n - 1) {
        throw std::invalid_argument("Diagonal size mismatch");
    }
    if (n <= 0) {
        throw std::invalid_argument("Matrix dimension must be positive");
    }
}

// Constructor for symmetric tridiagonal (off-diagonal elements are same)
TridiagonalMatrix::TridiagonalMatrix(const std::vector<double>& main,
                    const std::vector<double>& off_diag)
    : main_diag(main), n(main.size()) {
    
    if (static_cast<int>(off_diag.size()) != n - 1) {
        throw std::invalid_argument("Off-diagonal size mismatch");
    }
    
    lower_diag = off_diag;
    upper_diag = off_diag;
}

// Constructor for scalar tridiagonal (constant values)
TridiagonalMatrix::TridiagonalMatrix(int size, double lower_val, double main_val, double upper_val)
    : n(size) {
    
    if (size <= 0) {
        throw std::invalid_argument("Matrix size must be positive");
    }
    
    main_diag.assign(n, main_val);
    lower_diag.assign(n - 1, lower_val);
    upper_diag.assign(n - 1, upper_val);
}

TridiagonalMatrix::TridiagonalMatrix(const TridiagonalMatrix& other) 
    : lower_diag(other.lower_diag),
      main_diag(other.main_diag),
      upper_diag(other.upper_diag),
      n(other.n)
{
}

// Getters
int TridiagonalMatrix::size() const { return n; }
const std::vector<double>& TridiagonalMatrix::getLowerDiagonal() const { return lower_diag; }
const std::vector<double>& TridiagonalMatrix::getMainDiagonal() const { return main_diag; }
const std::vector<double>& TridiagonalMatrix::getUpperDiagonal() const { return upper_diag; }

// Element access (i, j)
double TridiagonalMatrix::operator()(int i, int j) const {
    if (i < 0 || i >= n || j < 0 || j >= n) {
        throw std::out_of_range("Matrix indices out of range");
    }
    
    if (i == j) {
        return main_diag[i];
    } else if (j == i + 1) {
        return upper_diag[i];
    } else if (j == i - 1) {
        return lower_diag[j];
    } else {
        return 0.0;  // Outside tridiagonal band
    }
}

// Set element (only within tridiagonal band)
void TridiagonalMatrix::set(int i, int j, double value) {
    if (i < 0 || i >= n || j < 0 || j >= n) {
        throw std::out_of_range("Matrix indices out of range");
    }
    
    if (i == j) {
        main_diag[i] = value;
    } else if (j == i + 1) {
        upper_diag[i] = value;
    } else if (j == i - 1) {
        lower_diag[j] = value;
    } else if (value != 0.0) {
        throw std::invalid_argument("Cannot set non-zero element outside tridiagonal band");
    }
}

// Solve Ax = b using Thomas algorithm
std::vector<double> TridiagonalMatrix::solve(const std::vector<double>& b) const {
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("Right-hand side vector size mismatch");
    }

    std::vector<double> c_prime(n - 1);  // Modified upper diagonal
    std::vector<double> d_prime(n);      // Modified RHS
    std::vector<double> x(n);            // Solution vector

    // Check for zero diagonal elements
    if (std::abs(main_diag[0]) < 1e-14) {
        throw std::runtime_error("Zero diagonal element encountered");
    }

    // Forward elimination
    c_prime[0] = upper_diag[0] / main_diag[0];
    d_prime[0] = b[0] / main_diag[0];

    for (int i = 1; i < n - 1; ++i) {
        double denom = main_diag[i] - lower_diag[i-1] * c_prime[i-1];
        
        if (std::abs(denom) < 1e-14) {
            throw std::runtime_error("Zero pivot encountered during elimination");
        }
        
        c_prime[i] = upper_diag[i] / denom;
        d_prime[i] = (b[i] - lower_diag[i-1] * d_prime[i-1]) / denom;
    }

    // Last row (no upper diagonal element)
    double denom = main_diag[n-1] - lower_diag[n-2] * c_prime[n-2];
    if (std::abs(denom) < 1e-14) {
        throw std::runtime_error("Zero pivot encountered in last row");
    }
    
    d_prime[n-1] = (b[n-1] - lower_diag[n-2] * d_prime[n-2]) / denom;

    // Back substitution
    x[n-1] = d_prime[n-1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i+1];
    }

    return x;
}

// Matrix-vector multiplication: y = Ax
std::vector<double> TridiagonalMatrix::multiply(const std::vector<double>& x) const {
    if (static_cast<int>(x.size()) != n) {
        throw std::invalid_argument("Vector size mismatch for multiplication");
    }

    std::vector<double> y(n, 0.0);

    // First row
    y[0] = main_diag[0] * x[0] + upper_diag[0] * x[1];

    // Middle rows
    for (int i = 1; i < n - 1; ++i) {
        y[i] = lower_diag[i-1] * x[i-1] + main_diag[i] * x[i] + upper_diag[i] * x[i+1];
    }

    // Last row
    y[n-1] = lower_diag[n-2] * x[n-2] + main_diag[n-1] * x[n-1];

    return y;
}

// Check if matrix is symmetric
bool TridiagonalMatrix::isSymmetric() const {
    for (int i = 0; i < n - 1; ++i) {
        if (std::abs(lower_diag[i] - upper_diag[i]) > 1e-14) {
            return false;
        }
    }
    return true;
}

bool TridiagonalMatrix::isPositiveDefinite() const {
    if (!isSymmetric()) {
        return false;
    }

    // For tridiagonal matrices, we can check this efficiently
    // by computing the Cholesky decomposition and seeing if it succeeds
    try {
        choleskyDecomposition();
        return true;
    } catch (const std::runtime_error&) {
        return false;
    }
}

std::pair<std::vector<double>, std::vector<double>> TridiagonalMatrix::choleskyDecomposition() const {
    if (!isSymmetric()) {
        throw std::runtime_error("Cholesky decomposition requires symmetric matrix");
    }

    std::vector<double> L_diag(n);      // Diagonal elements of L
    std::vector<double> L_lower(n-1);   // Lower diagonal elements of L

    // First element
    if (main_diag[0] <= 0) {
        throw std::runtime_error("Matrix is not positive definite");
    }
    L_diag[0] = std::sqrt(main_diag[0]);

    // Compute L diagonal and lower diagonal elements
    for (int i = 1; i < n; ++i) {
        // Lower diagonal element: L[i][i-1] = A[i][i-1] / L[i-1][i-1]
        L_lower[i-1] = lower_diag[i-1] / L_diag[i-1];
        
        // Diagonal element: L[i][i] = sqrt(A[i][i] - L[i][i-1]^2)
        double temp = main_diag[i] - L_lower[i-1] * L_lower[i-1];
        if (temp <= 0) {
            std::cout << "here" << std::endl;
            this->print();
            throw std::runtime_error("Matrix is not positive definite");
        }
        L_diag[i] = std::sqrt(temp);
    }

    return std::make_pair(L_diag, L_lower);
}

// Solve using Cholesky decomposition: A = L * L^T, solve L * y = b, then L^T * x = y
std::vector<double> TridiagonalMatrix::solveCholesky(const std::vector<double>& b) const {
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("Right-hand side vector size mismatch");
    }

    // Get Cholesky decomposition
    auto [L_diag, L_lower] = choleskyDecomposition();

    std::vector<double> y(n);  // Intermediate solution
    std::vector<double> x(n);  // Final solution

    // Forward substitution: L * y = b
    y[0] = b[0] / L_diag[0];
    for (int i = 1; i < n; ++i) {
        y[i] = (b[i] - L_lower[i-1] * y[i-1]) / L_diag[i];
    }

    // Back substitution: L^T * x = y
    x[n-1] = y[n-1] / L_diag[n-1];
    for (int i = n-2; i >= 0; --i) {
        x[i] = (y[i] - L_lower[i] * x[i+1]) / L_diag[i];
    }

    return x;
}



// Print matrix (for debugging)
void TridiagonalMatrix::print(std::ostream& os) const {
    os << std::fixed << std::setprecision(3);
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            os << std::setw(8) << (*this)(i, j) << " ";
        }
        os << std::endl;
    }
}

std::vector<double> solveLowerBidiagonal(
    const std::pair<std::vector<double>, std::vector<double>>& L,
    const std::vector<double>& b) {
    
    const auto& diagonal = L.first;
    const auto& subdiagonal = L.second;
    int n = diagonal.size();
    
    if (static_cast<int>(subdiagonal.size()) != n - 1) {
        throw std::invalid_argument("Subdiagonal size must be n-1");
    }
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("Right-hand side vector size mismatch");
    }

    std::vector<double> x(n);

    // Forward substitution for lower bidiagonal system
    // L[0,0] * x[0] = b[0]
    if (std::abs(diagonal[0]) < 1e-14) {
        throw std::runtime_error("Zero diagonal element encountered");
    }
    x[0] = b[0] / diagonal[0];

    // For i > 0: L[i,i] * x[i] + L[i,i-1] * x[i-1] = b[i]
    // Therefore: x[i] = (b[i] - L[i,i-1] * x[i-1]) / L[i,i]
    for (int i = 1; i < n; ++i) {
        if (std::abs(diagonal[i]) < 1e-14) {
            throw std::runtime_error("Zero diagonal element encountered");
        }
        x[i] = (b[i] - subdiagonal[i-1] * x[i-1]) / diagonal[i];
    }

    return x;
}

std::vector<double> solveUpperBidiagonal(
    const std::pair<std::vector<double>, std::vector<double>>& U,
    const std::vector<double>& b) {
    
    const auto& diagonal = U.first;
    const auto& superdiagonal = U.second;
    int n = diagonal.size();
    
    if (static_cast<int>(superdiagonal.size()) != n - 1) {
        throw std::invalid_argument("Superdiagonal size must be n-1");
    }
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("Right-hand side vector size mismatch");
    }

    std::vector<double> x(n);

    // Back substitution for upper bidiagonal system
    // U[n-1,n-1] * x[n-1] = b[n-1]
    if (std::abs(diagonal[n-1]) < 1e-14) {
        throw std::runtime_error("Zero diagonal element encountered");
    }
    x[n-1] = b[n-1] / diagonal[n-1];

    // For i < n-1: U[i,i] * x[i] + U[i,i+1] * x[i+1] = b[i]
    // Therefore: x[i] = (b[i] - U[i,i+1] * x[i+1]) / U[i,i]
    for (int i = n-2; i >= 0; --i) {
        if (std::abs(diagonal[i]) < 1e-14) {
            throw std::runtime_error("Zero diagonal element encountered");
        }
        x[i] = (b[i] - superdiagonal[i] * x[i+1]) / diagonal[i];
    }

    return x;
}

void demonstrateTridiagonalMatrix() {
    std::cout << "Tridiagonal Matrix Class Demonstration" << std::endl;
    std::cout << "=====================================" << std::endl;

    try {
        // Example 1: General tridiagonal matrix
        std::vector<double> lower = {1.0, 1.0, 1.0};  // Lower diagonal
        std::vector<double> main = {2.0, 2.0, 2.0, 2.0};  // Main diagonal
        std::vector<double> upper = {1.0, 1.0, 1.0}; // Upper diagonal
        
        TridiagonalMatrix A(lower, main, upper);
        
        std::cout << "Matrix A:" << std::endl;
        A.print();
        std::cout << std::endl;

        // Test system: Ax = b
        std::vector<double> b = {1.0, 4.0, 7.0, 10.0};
        std::cout << "Right-hand side b: ";
        for (double val : b) std::cout << val << " ";
        std::cout << std::endl << std::endl;

        // Solve the system
        std::vector<double> x = A.solve(b);
        
        std::cout << "Solution x:" << std::endl;
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            std::cout << "x[" << i << "] = " << std::fixed << std::setprecision(6) 
                      << x[i] << std::endl;
        }
        std::cout << std::endl;

        // Verify: compute Ax and compare with b
        std::vector<double> Ax = A.multiply(x);
        std::cout << "Verification (Ax should equal b):" << std::endl;
        for (int i = 0; i < static_cast<int>(Ax.size()); ++i) {
            std::cout << "Ax[" << i << "] = " << std::fixed << std::setprecision(6) 
                      << Ax[i] << ", b[" << i << "] = " << b[i] 
                      << ", error = " << std::abs(Ax[i] - b[i]) << std::endl;
        }
        std::cout << std::endl;

        // Example 2: Symmetric tridiagonal matrix
        std::cout << "Example 2: Symmetric tridiagonal matrix" << std::endl;
        std::vector<double> main_sym = {4.0, 4.0, 4.0};
        std::vector<double> off_diag = {-1.0, -1.0};
        
        TridiagonalMatrix B(main_sym, off_diag);
        std::cout << "Matrix B (symmetric):" << std::endl;
        B.print();
        std::cout << "Is symmetric: " << (B.isSymmetric() ? "Yes" : "No") << std::endl;
        std::cout << "Is positive definite: " << (B.isPositiveDefinite() ? "Yes" : "No") << std::endl;
        std::cout << std::endl;

        // Test Cholesky decomposition
        if (B.isPositiveDefinite()) {
            std::cout << "Testing Cholesky decomposition:" << std::endl;
            auto [L_diag, L_lower] = B.choleskyDecomposition();
            
            std::cout << "Cholesky factor L diagonal: ";
            for (double val : L_diag) std::cout << std::fixed << std::setprecision(4) << val << " ";
            std::cout << std::endl;
            
            std::cout << "Cholesky factor L lower: ";
            for (double val : L_lower) std::cout << std::fixed << std::setprecision(4) << val << " ";
            std::cout << std::endl;

            // Test solve with Cholesky
            std::vector<double> b_test = {1.0, 2.0, 3.0};
            std::vector<double> x_thomas = B.solve(b_test);
            std::vector<double> x_cholesky = B.solveCholesky(b_test);
            
            std::cout << "\nComparing Thomas vs Cholesky solutions:" << std::endl;
            for (int i = 0; i < static_cast<int>(x_thomas.size()); ++i) {
                std::cout << "x[" << i << "]: Thomas = " << std::fixed << std::setprecision(6) 
                          << x_thomas[i] << ", Cholesky = " << x_cholesky[i] 
                          << ", diff = " << std::abs(x_thomas[i] - x_cholesky[i]) << std::endl;
            }
        }

        // Example 3: Test bidiagonal solvers
        std::cout << "\nExample 3: Testing bidiagonal solvers" << std::endl;
        
        // Test lower bidiagonal system
        std::vector<double> L_diag = {2.0, 3.0, 1.0};      // diagonal
        std::vector<double> L_sub = {1.0, 0.5};             // subdiagonal
        std::vector<double> b_lower = {4.0, 7.0, 2.0};
        
        std::cout << "Lower bidiagonal system L:" << std::endl;
        std::cout << "[ 2.0  0.0  0.0 ]" << std::endl;
        std::cout << "[ 1.0  3.0  0.0 ]" << std::endl;
        std::cout << "[ 0.0  0.5  1.0 ]" << std::endl;
        std::cout << "b = [4.0, 7.0, 2.0]" << std::endl;
        
        auto L_matrix = std::make_pair(L_diag, L_sub);
        std::vector<double> x_lower = solveLowerBidiagonal(L_matrix, b_lower);
        
        std::cout << "Solution:" << std::endl;
        for (int i = 0; i < static_cast<int>(x_lower.size()); ++i) {
            std::cout << "x[" << i << "] = " << std::fixed << std::setprecision(6) << x_lower[i] << std::endl;
        }
        
        // Test upper bidiagonal system
        std::vector<double> U_diag = {2.0, 3.0, 1.0};      // diagonal
        std::vector<double> U_super = {1.0, 0.5};           // superdiagonal
        std::vector<double> b_upper = {5.0, 8.5, 2.0};
        
        std::cout << "\nUpper bidiagonal system U:" << std::endl;
        std::cout << "[ 2.0  1.0  0.0 ]" << std::endl;
        std::cout << "[ 0.0  3.0  0.5 ]" << std::endl;
        std::cout << "[ 0.0  0.0  1.0 ]" << std::endl;
        std::cout << "b = [5.0, 8.5, 2.0]" << std::endl;
        
        auto U_matrix = std::make_pair(U_diag, U_super);
        std::vector<double> x_upper = solveUpperBidiagonal(U_matrix, b_upper);
        
        std::cout << "Solution:" << std::endl;
        for (int i = 0; i < static_cast<int>(x_upper.size()); ++i) {
            std::cout << "x[" << i << "] = " << std::fixed << std::setprecision(6) << x_upper[i] << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}