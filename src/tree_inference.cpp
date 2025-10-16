#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <string>
#include <vector>
#include "sim.h"
#include "random.h"
#include "tree_inference.h"


State constructInitialState(
    const PerfectPhylo& perfectPhylo,
    double ntInit // Initial value of nt, treated as constant in initial state
) {

    State state;
    state.n = perfectPhylo.n;
    state.m = perfectPhylo.m;
    state.nt.resize(state.m+state.n-1, ntInit);
    state.g.resize(state.m+state.n-1, -1);
    state.produced.resize(state.m+state.n-1, -1);
    state.ng.resize(state.m+state.n-1, -1);
    state.k.resize(state.m+state.n-1, -1);
    state.c.resize(state.m+state.n-1, -1);
    state.rooted = perfectPhylo.rooted;

    // What index of the evolutionary sequence are we on?
    int currIdx = 0;
    int k = state.n;

    // Loop backwards over nodes in perfect phylo to obtain post-order traversal
    for(int g = state.m; g >= 0; g--) {

        // Number of extant lineages when g is called
        int ng = perfectPhylo.w[g] + perfectPhylo.neighbors[g].size() - 1;

        // If g is the root, this increases by 1
        if(g == 0) {
            ng++;
        }

        // Coalescences of genotype g
        while(ng > 1) {
            state.g[currIdx] = g;
            state.ng[currIdx] = ng;
            state.k[currIdx] = k;
            state.c[currIdx] = k * (k-1) / 2;

            k--; // Total number of extant lineages decreases
            ng--;
            currIdx++;
        }

        // If g is not the root, do a mutation
        if(g != 0) {

            if(perfectPhylo.neighbors[g].empty()) {
                throw std::runtime_error("Expected at least one neighbor for non-root genotype");
            }
            state.g[currIdx] = g;
            state.ng[currIdx] = ng;
            state.k[currIdx] = k;
            state.c[currIdx] = k * (k-1) / 2;
            state.produced[currIdx] = perfectPhylo.neighbors[g][0]; // Ancestral genotype is always the first neighbor!

            ng--;
            // Nothing happens to k
            currIdx++;
        }
    }

    if(currIdx != state.m+state.n-1) {
        std::cout << "currIdx: " << currIdx << std::endl;
        throw std::runtime_error("Called wrong number of times");
    }

    return state;
}

void printState(const State& state) {
    std::cout << "=== State ===" << std::endl;
    std::cout << "Number of leaves (n): " << state.n << std::endl;
    std::cout << "Total mutations (m): " << state.m << std::endl;
    std::cout << "Rooted: " << (state.rooted ? "yes" : "no") << std::endl;
    std::cout << std::endl;
    
    std::cout << "Vector sizes:" << std::endl;
    std::cout << "  nt.size(): " << state.nt.size() << std::endl;
    std::cout << "  g.size(): " << state.g.size() << std::endl;
    std::cout << "  produced.size(): " << state.produced.size() << std::endl;
    std::cout << "  ng.size(): " << state.ng.size() << std::endl;
    std::cout << "  k.size(): " << state.k.size() << std::endl;
    std::cout << "  c.size(): " << state.c.size() << std::endl;
    std::cout << std::endl;
    
    std::cout << "Effective population sizes (nt):" << std::endl;
    for (size_t i = 0; i < state.nt.size(); i++) {
        std::cout << "  [" << i << "]: " << state.nt[i] << std::endl;
    }
    std::cout << std::endl;
    
    // Use the minimum size to avoid out-of-bounds access
    size_t numEvents = state.g.size();
    
    std::cout << "Event information:" << std::endl;
    for (size_t i = 0; i < numEvents; i++) {
        std::cout << "  Event " << i << ":" << std::endl;
        std::cout << "    Genotype: " << state.g[i] << std::endl;
        std::cout << "    Produced genotype: " << state.produced[i];
        if (state.produced[i] == -1) {
            std::cout << " (coalescence)";
        } else {
            std::cout << " (mutation)";
        }
        std::cout << std::endl;
        std::cout << "    Extant lineages of genotype g prior: " << state.ng[i] << std::endl;
        std::cout << "    Total extant lineages prior (k): " << state.k[i] << std::endl;
        std::cout << "    k choose 2 (c): " << state.c[i] << std::endl;
        std::cout << std::endl;
    }
}