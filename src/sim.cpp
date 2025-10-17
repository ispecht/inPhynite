#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <string>
#include "sim.h"
#include "random.h"
#include "helpers.h"

Tree simulateTree(
    int n,
    double ntInit,
    double ntLogSlope,
    std::mt19937& rng
) {

    Tree tree;

    tree.n = n;
    tree.m = 0;
    tree.nt.clear();
    tree.h.resize(2*n - 1, -1);
    tree.children.resize(2*n - 1, {});
    tree.t.resize(2*n - 1, 0.0);
    tree.tMut.resize(2*n - 1, {});

    // Branches that can still coalesce
    std::vector<int> active(n);
    std::iota(active.begin(), active.end(), 0);

    double tCurr = 0.0;
    double ntCurr = ntInit;
    int nActive = n;

    for(int i = n; i < 2*n - 1; i++) {

        tree.nt.push_back(ntCurr);

        double coalRate = (double) (nActive * (nActive - 1) / 2);
        coalRate /= ntCurr;

        // Update current time
        double deltaT = rExp(coalRate, rng);
        tCurr += deltaT;

        // Update ntCurr
        ntCurr = std::exp(
            std::log(ntCurr) + (ntLogSlope * deltaT)
        );

        nActive--;

        // Time of node i
        tree.t[i] = tCurr;

        // Sample the children
        std::pair<int, int> children = rPair(active, rng);
        tree.children[i].push_back(children.first);
        tree.children[i].push_back(children.second);

        for(int j : tree.children[i]) {
            // Update active set
            removeValue(active, j);

            // Ancestor of j
            tree.h[j] = i;

            // Branch length
            double len = tree.t[i] - tree.t[j];

            if(len <= 0.0) {
                throw std::runtime_error("Non-positive branch length");
            }

            // Number of mutations along branch j
            int nMut = rPois(len, rng);
            tree.m += nMut;

            // Times of mutations along said branch
            tree.tMut[j] = rUnifSorted(tree.t[j], tree.t[i], nMut, rng);
        }

        // Update active set
        active.push_back(i);

    }

    return tree;
}

void extractPerfectPhyloRec(
    const Tree& tree,
    int treeNode,
    int ancGeno, // genotype of the ancestral node
    int& nextAvail, // Next available genotype to be filled out
    PerfectPhylo& soFar
) {

    // The last node in the perfect phylo is the current genotype
    for(int j = 0; j < (int) tree.tMut[treeNode].size(); j++) {
        soFar.neighbors[ancGeno].push_back(nextAvail);
        soFar.neighbors[nextAvail].push_back(ancGeno);
        
        ancGeno = nextAvail;
        nextAvail++;
    }

    std::vector<int> children = tree.children[treeNode];

    // ancGeno is now the genotype of treeNode
    if(children.size() == 0) {
        // treeNode is a leaf
        soFar.w[ancGeno]++;
        soFar.leaves[ancGeno].push_back(treeNode);
        return;
    }

    for(int j : children) {
        extractPerfectPhyloRec(tree, j, ancGeno, nextAvail, soFar);
    }
}


PerfectPhylo extractPerfectPhylo(
    const Tree& tree,
    bool rooted
) {

    PerfectPhylo perfectPhylo;
    perfectPhylo.n = tree.n;
    perfectPhylo.m = tree.m;
    perfectPhylo.neighbors.resize(tree.m + 1, {});
    perfectPhylo.leaves.resize(tree.m + 1, {});
    perfectPhylo.w.resize(tree.m + 1, 0);
    perfectPhylo.rooted = rooted;


    int nextAvail = 1;
    extractPerfectPhyloRec(tree, 2*tree.n - 2, 0, nextAvail, perfectPhylo);

    return perfectPhylo;

}

void printPerfectPhylo(const PerfectPhylo& phylo) {
    std::cout << "=== Perfect Phylogeny ===" << std::endl;
    std::cout << "Number of leaves (n): " << phylo.n << std::endl;
    std::cout << "Number of mutations (m): " << phylo.m << std::endl;
    std::cout << "Rooted?: " << phylo.rooted << std::endl;
    std::cout << std::endl;
    
    std::cout << "Node information:" << std::endl;
    for (size_t i = 0; i < phylo.w.size(); i++) {
        std::cout << "  Node " << i << ":" << std::endl;
        std::cout << "    Weight: " << phylo.w[i] << std::endl;
        
        std::cout << "    Neighbors: ";
        if (phylo.neighbors[i].empty()) {
            std::cout << "none";
        } else {
            for (size_t j = 0; j < phylo.neighbors[i].size(); j++) {
                std::cout << phylo.neighbors[i][j];
                if (j < phylo.neighbors[i].size() - 1) std::cout << ", ";
            }
        }

        std::cout << std::endl;

        std::cout << "    Leaves: ";
        if (phylo.leaves[i].empty()) {
            std::cout << "none";
        } else {
            for (size_t j = 0; j < phylo.leaves[i].size(); j++) {
                std::cout << phylo.leaves[i][j];
                if (j < phylo.leaves[i].size() - 1) std::cout << ", ";
            }
        }

        std::cout << std::endl;
        std::cout << std::endl;
    }
}













void printTree(const Tree& tree) {
    std::cout << "=== Tree Structure ===" << std::endl;
    std::cout << "Number of leaves (n): " << tree.n << std::endl;
    std::cout << "Total mutations (m): " << tree.m << std::endl;
    std::cout << std::endl;
    
    std::cout << "Effective population sizes (nt):" << std::endl;
    for (size_t i = 0; i < tree.nt.size(); i++) {
        std::cout << "  [" << i << "]: " << tree.nt[i] << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Node information:" << std::endl;
    for (size_t i = 0; i < tree.h.size(); i++) {
        std::cout << "  Node " << i << ":" << std::endl;
        std::cout << "    Parent: " << tree.h[i] << std::endl;
        std::cout << "    Time: " << tree.t[i] << std::endl;
        
        std::cout << "    Children: ";
        if (tree.children[i].empty()) {
            std::cout << "none (leaf)";
        } else {
            for (size_t j = 0; j < tree.children[i].size(); j++) {
                std::cout << tree.children[i][j];
                if (j < tree.children[i].size() - 1) std::cout << ", ";
            }
        }
        std::cout << std::endl;
        
        std::cout << "    Mutations: ";
        if (tree.tMut[i].empty()) {
            std::cout << "none";
        } else {
            std::cout << tree.tMut[i].size() << " mutation(s) at times: ";
            for (size_t j = 0; j < tree.tMut[i].size(); j++) {
                std::cout << tree.tMut[i][j];
                if (j < tree.tMut[i].size() - 1) std::cout << ", ";
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
}