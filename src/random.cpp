#include "random.h"
#include <algorithm>
#include <cmath>

// Poisson random variable
int rPois(double lambda, std::mt19937& rng) {
    std::poisson_distribution<int> dist(lambda);
    return dist(rng);
}

// Exponential random variable
double rExp(double lambda, std::mt19937& rng) {
    std::exponential_distribution<double> dist(lambda);
    return dist(rng);
}

// Uniform random variable
double rUnif(double a, double b, std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

std::vector<double> rUnifSorted(double a, double b, int n, std::mt19937& rng) {
    // Input validation
    if (a > b) {
        throw std::invalid_argument("Parameter 'a' must be less than 'b'");
    }
    if (n < 0) {
        throw std::invalid_argument("Parameter 'n' must be nonnegative");
    }

    if(n == 0){
        return {};
    }

    // Set up random number generation
    std::uniform_real_distribution<double> dist(a, b);
    
    // Generate n random numbers
    std::vector<double> result;
    result.reserve(n);
    
    for (int i = 0; i < n; ++i) {
        result.push_back(dist(rng));
    }
    
    // Sort the vector
    std::sort(result.begin(), result.end());
    
    return result;
}

// Normal random variable
double rNorm(double mu, double sigma, std::mt19937& rng) {
    std::normal_distribution<double> dist(mu, sigma);
    return dist(rng);
}

// Sample two elements without replacement
std::pair<int, int> rPair(const std::vector<int>& elems, std::mt19937& rng) {
    if (elems.size() < 2) {
        throw std::invalid_argument("rPair requires at least 2 elements");
    }
    
    std::uniform_int_distribution<size_t> dist1(0, elems.size() - 1);
    size_t idx1 = dist1(rng);
    
    std::uniform_int_distribution<size_t> dist2(0, elems.size() - 2);
    size_t idx2 = dist2(rng);
    
    // Adjust idx2 if it's >= idx1
    if (idx2 >= idx1) {
        idx2++;
    }
    
    return std::make_pair(elems[idx1], elems[idx2]);
}