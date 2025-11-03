#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <random>
#include "sim.h"

// Structure for storing phylogenetic tree annotated with mutations
struct State {
    int n; // n is always number of leaves
    int m; // total number of mutations
    std::vector<double> nt; // Effective population size on each intercoalescent interval
    std::vector<int> g; // Genotype in which each mutation/coalescent occurs
    std::vector<int> produced; // Genotype of new lineage produced via mutation. -1 if coalescence
    std::vector<int> ng; // Number of extant lineages of genotype g immediately PRIOR to its coalescence/mutation
    std::vector<int> k; // Number of TOTAL extant lineages immediately PRIOR to ith event
    std::vector<int> c; // k choose 2; stored for the sake of efficiency
    bool rooted; // Is the root genotype fixed?
};

State constructInitialState(
    const PerfectPhylo& perfectPhylo,
    double ntInit // Initial value of nt, treated as constant in initial state
);

void updateTree(
    State& state,
    int i,
    std::mt19937& rng
);

void printState(const State& state);
