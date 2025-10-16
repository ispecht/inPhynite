#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <random>

// Poisson
int rPois(double lambda, std::mt19937& rng);

// Expo
double rExp(double lambda, std::mt19937& rng);

// Uniform
double rUnif(double a, double b, std::mt19937& rng);

// Sorted n uniforms low to high
std::vector<double> rUnifSorted(double a, double b, int n, std::mt19937& rng);

// Normal
double rNorm(double mu, double sigma, std::mt19937& rng);

// Sample two elements from a vector of elements without replacement
std::pair<int, int> rPair(const std::vector<int>& elems, std::mt19937& rng);
