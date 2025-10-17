#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <random>

// Structure for storing phylogenetic tree annotated with mutations
struct Tree {
    int n; // n is always number of leaves
    int m; // total number of mutations
    std::vector<double> nt; // Effective population size on each intercoalescent interval
    std::vector<int> h; // parent
    std::vector<std::vector<int>> children; // children
    std::vector<double> t; // time of each node
    std::vector<std::vector<double>> tMut; // time of mutations along each branch
};

// Augmented perfect phylogeny
struct PerfectPhylo {
    int n; // total number of leaves
    int m; // total number of mutations; equals number of nodes in the perfect phylo minus one
    std::vector<int> w; // Weight of each node
    std::vector<std::vector<int>> neighbors; // Neighbors of each node
    std::vector<std::vector<int>> leaves; // Leaves corresponding to each genotype
    bool rooted; // If rooted, the root genotype is always index 0
};

Tree simulateTree(
    int n,
    double ntInit,
    double ntLogSlope,
    std::mt19937& rng
);

void printTree(const Tree& tree);

PerfectPhylo extractPerfectPhylo(
    const Tree& tree,
    bool rooted
);

void printPerfectPhylo(const PerfectPhylo& phylo);