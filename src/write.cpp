#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <string>
#include "sim.h"
#include "random.h"
#include "write.h"
#include "sim.h"
#include "tree_inference.h"
#include "helpers.h"


Tree stateToTree(
    const State& state,
    PerfectPhylo perfectPhylo, // passing by value because we will update this
    std::mt19937& rng
) {

    Tree tree;

    tree.n = state.n;
    tree.m = state.m;
    tree.nt = state.nt;
    tree.h.resize(2*state.n - 1, -1);
    tree.children.resize(2*state.n - 1, {});
    tree.t.resize(2*state.n - 1, 0.0);
    tree.tMut.resize(2*state.n - 1, {}); // Ignoring this

    double tCurr = 0.0;

    int nextTreeNodeIdx = state.n;

    for(int i = 0; i < state.m + state.n - 1; i++) {

        double rate = (double) state.k[i] + (((double) state.c[i]) / state.nt[state.n - state.k[i]]);

        // Update current time
        double deltaT = rExp(rate, rng);
        tCurr += deltaT;

        int g = state.g[i];

        // Check agreement with ng
        if((int) perfectPhylo.leaves[g].size() != state.ng[i]) {
            throw std::runtime_error("ng disagreement");
        }

        // If mutation, update perfect phylo only
        if(state.produced[i] != -1) {

            // Should be exactly one extant lineage of genotype g
            if(perfectPhylo.leaves[g].size() != 1) {
                throw std::runtime_error("Should be exactly one extant lineage of mutating genotype");
            }

            int treeNodeIdx = perfectPhylo.leaves[g][0];
            perfectPhylo.leaves[state.produced[i]].push_back(treeNodeIdx);

        }else{

            // Should be 2+ extant lineages of genotype g
            if(perfectPhylo.leaves[g].size() <= 1) {
                throw std::runtime_error("Should be 2+ extant lineages of mutating genotype");
            }

            // Time of node nextTreeNodeIdx
            tree.t[nextTreeNodeIdx] = tCurr;

            // Sample the children
            std::pair<int, int> children = rPair(perfectPhylo.leaves[g], rng);
            tree.children[nextTreeNodeIdx].push_back(children.first);
            tree.children[nextTreeNodeIdx].push_back(children.second);

            for(int j : tree.children[nextTreeNodeIdx]) {
                // Update active set
                removeValue(perfectPhylo.leaves[g], j);

                // Ancestor of j
                tree.h[j] = nextTreeNodeIdx;

            }

            perfectPhylo.leaves[g].push_back(nextTreeNodeIdx);


            nextTreeNodeIdx++;

        }
    }

    if(nextTreeNodeIdx != 2*state.n - 1) {
        throw std::runtime_error("wrong end node index");
    }

    return tree;

}


std::string treeToNewickRec(const Tree& tree, int node) {
    // If leaf node
    if (tree.children[node].empty()) {
        return std::to_string(node);
    }
    
    // Internal node - process children
    std::string result = "(";
    for (size_t i = 0; i < tree.children[node].size(); i++) {
        int child = tree.children[node][i];
        result += treeToNewickRec(tree, child);
        
        // Branch length is difference in times
        double branchLength = tree.t[node] - tree.t[child];
        result += ":" + std::to_string(branchLength);
        
        if (i < tree.children[node].size() - 1) {
            result += ",";
        }
    }
    result += ")";
    
    // Add node label for internal nodes
    result += std::to_string(node);
    
    return result;
}

std::string treeToNewick(const Tree& tree) {
    // Find root (node with no parent, h[i] == -1)
    int root = -1;
    for (size_t i = 0; i < tree.h.size(); i++) {
        if (tree.h[i] == -1) {
            root = i;
            break;
        }
    }
    
    if (root == -1) {
        throw std::runtime_error("No root found in tree");
    }
    
    return treeToNewickRec(tree, root) + ";";
}


std::vector<std::string> perfectPhyloToSequences(const PerfectPhylo& perfectPhylo, int L) {

    // Idea: for each node on perfect phylo in ancestral history of sequence, we get a mutation at the index of that node

    if(L < (int) perfectPhylo.w.size()) {
        throw std::runtime_error("L too small");
    }

    std::string rootSeq(L, 'A');
    std::vector<std::string> sequences(perfectPhylo.n, rootSeq);

    for(int g = 0; g < (int) perfectPhylo.w.size(); g++) {
        // Loop over hosts of that genotype
        for(int j : perfectPhylo.leaves[g]) {
            // Current genotype
            int gCurr = g;
            while(gCurr != 0) {
                sequences[j][gCurr] = 'C';
                gCurr = perfectPhylo.neighbors[gCurr][0];
            }
        }
    }

    return sequences;

}

