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
    state.nt.resize(state.n-1, ntInit);
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

void updateTree(
    State& state,
    int i,
    std::mt19937& rng
) {

    // Smaller index of neighbor swap
    //int i = (int) rUnif(0.0, state.m + state.n - 2, rng);

    // If genotypes are the same, nothing to do here!
    if(state.g[i] == state.g[i+1]) {
        return;
    }

    // If two coalescences, auto accept
    if(
        (state.produced[i] == -1) &&
        (state.produced[i+1] == -1)
    ) {
        std::swap(state.g[i], state.g[i+1]);
        std::swap(state.ng[i], state.ng[i+1]);
        return;
    }

    // If two mutations, and the first doesn't produce the second genotype, auto accept
    if(
        (state.produced[i] != -1) &&
        (state.produced[i+1] != -1)
    ) {
        if(state.produced[i] != state.g[i+1]) {
            std::swap(state.g[i], state.g[i+1]);
            std::swap(state.produced[i], state.produced[i+1]);
            std::swap(state.ng[i], state.ng[i+1]);
            return;
        }
    }

    // Check if we're affecting the root
    bool affectsRoot = false;
    if(i == state.m + state.n - 3) {
        affectsRoot = true;
    }

    // Check if the move is valid; if not, auto reject

    if(affectsRoot) {
        // If the tree is rooted, auto reject
        if(state.rooted) {
            return;
        }

        // Otherwise, the swap is OK if it's exchanging a mutation for a coalescence
        // First, check last index is always a coalescence
        if(state.produced[i+1] != -1) {
            throw std::runtime_error("Last event should always be a coalescence");
        }

        // If both are coalescences, auto reject, move makes no impact on tree
        if(state.produced[i] == -1) {
            return;
        }

        // If mutation then coalescence and tree is unrooted, can accept the move, and no impact on likelihood
        // Be careful to update all elements of state correctly

        std::swap(state.g[i], state.g[i+1]);
        // ng doesn't change: the second to last genotype must have 1 extant lineage, the last must have 2
        state.produced[i] = state.g[i+1]; // Now the updated root genotype is state.g[i+1]
        // state.produced[i+1] is already -1

        // Likelihood doesn't change; we're done
        return;
    }

    // The root is not being affected if we got here

    // The scenario in which we must auto-reject the move is:
    // ith event is a mutation
    // ith event produces (i+1)st genotype
    // (i+1)st event has 2 or fewer extant lineages
    if(state.produced[i] == state.g[i+1]) {
        if(state.ng[i+1] <= 2) {
            // Auto reject
            return;
        }
    }

    // If we got here, we should be in the regime of coal, mut or mut, coal
    // Remember that the genotypes of the two events are necessarily different

    if(
        (state.produced[i] == -1) &&
        (state.produced[i+1] != -1)
    ) {
        // Current state: coalescence then mutation

        // c[i], k[i] stay the same
        // c[i+1], k[i+1] change
        int curr_k = state.k[i+1];
        int curr_c = state.c[i+1];

        // Old effective population size on interval we care about
        double curr_nt = state.nt[state.n - curr_k];

        int prop_k = curr_k + 1;
        int prop_c = prop_k * (prop_k - 1) / 2;

        // New effective population size on interval we care about
        double prop_nt = state.nt[state.n - prop_k];

        // Current number of lineages of coalescing genotype
        int curr_ng = state.ng[i];

        // posterior log ratio
        double logRatio = 0.0;

        logRatio -= std::log((double)prop_k + ((double)prop_c / prop_nt));
        logRatio += std::log((double)curr_k + ((double)curr_c / curr_nt));


        // Check if we pick up an extra coalescing lineage in new state
        int prop_ng = curr_ng;
        if(state.produced[i+1] == state.g[i]) {
            prop_ng++;

            // Factor this into calculation
            // Posterior ratio:
            // (prop_ng choose 2) / (curr_ng choose 2) 
            // = ((prop_ng)(prop_ng - 1)) / ((curr_ng)(curr_ng-1))
            // = prop_ng / (curr_ng - 1)
            logRatio += std::log((double)prop_ng);
            logRatio -= std::log((double)curr_ng - 1.0);
        }

        // MH accept/reject criterion
        double u = rUnif(0.0, 1.0, rng);

        if(std::log(u) < logRatio) {

            std::swap(state.g[i], state.g[i+1]);
            std::swap(state.produced[i], state.produced[i+1]);
            std::swap(state.ng[i], state.ng[i+1]);
            if(state.produced[i] == state.g[i+1]) {
                state.ng[i+1]++;
            }
            state.k[i+1] = prop_k;
            state.c[i+1] = prop_c;
        }

    }else if(
        (state.produced[i] != -1) &&
        (state.produced[i+1] == -1)
    ) {
        // Current state: mutation then coalescence

        // c[i], k[i] stay the same
        // c[i+1], k[i+1] change
        int curr_k = state.k[i+1];
        int curr_c = state.c[i+1];

        // Old effective population size on interval we care about
        double curr_nt = state.nt[state.n - curr_k];

        int prop_k = curr_k - 1; // NOW: this decreases by 1
        int prop_c = prop_k * (prop_k - 1) / 2;

        // New effective population size on interval we care about
        double prop_nt = state.nt[state.n - prop_k];

        // Current number of lineages of coalescing genotype
        int curr_ng = state.ng[i+1]; // NOW: i+1 not i

        if(state.produced[i] == state.g[i+1]) {
            if(curr_ng <= 2) {
                throw std::runtime_error("Should be at least 3");
            }
        }
        

        // posterior log ratio
        double logRatio = 0.0;

        logRatio -= std::log((double)prop_k + ((double)prop_c / prop_nt));
        logRatio += std::log((double)curr_k + ((double)curr_c / curr_nt));


        // Check if we pick up an extra coalescing lineage in CURRENT state
        int prop_ng = curr_ng;
        if(state.produced[i] == state.g[i+1]) {
            prop_ng--;

            // Factor this into calculation
            // Posterior ratio:
            // (prop_ng choose 2) / (curr_ng choose 2) 
            // = ((prop_ng)(prop_ng - 1)) / ((curr_ng)(curr_ng-1))
            // = (prop_ng-1) / curr_ng
            logRatio += std::log((double)prop_ng - 1.0);
            logRatio -= std::log((double)curr_ng);
        }

        // MH accept/reject criterion
        double u = rUnif(0.0, 1.0, rng);

        if(std::log(u) < logRatio) {

            std::swap(state.g[i], state.g[i+1]);
            std::swap(state.produced[i], state.produced[i+1]);
            std::swap(state.ng[i], state.ng[i+1]);
            if(state.produced[i+1] == state.g[i]) { // NOW: REVERSED
                state.ng[i]--;
            }
            state.k[i+1] = prop_k;
            state.c[i+1] = prop_c;
        }
    }else{
        std::cout << "Debug:" << std::endl;
        std::cout << state.produced[i] << std::endl;
        std::cout << state.produced[i+1] << std::endl;
        throw std::runtime_error("Should be in coal, mut or mut, coal regime");
    }
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