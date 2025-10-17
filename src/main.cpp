#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include "sim.h"
#include "tree_inference.h"
#include "nt_inference.h"
#include "write.h"


int main() {

    double NT_INIT = 32.0;
    double NT_LOG_SLOPE = 0.0;
    int N = 40;
    int N_GLOBAL_ITERS = 10000;
    int N_TREE_ITERS = 10000;
    int SAMPLE_EVERY = 10;
    bool INFER_NT = true;
    std::string OUT_DIR = "/Users/ispecht/Desktop/demo_output";
    const int SEED = 123;



    std::mt19937 rng;
    rng.seed(SEED);

    Tree tree = simulateTree(
        N,
        NT_INIT,
        NT_LOG_SLOPE,
        rng
    );

    std::cout << "Number of mutations on tree: " << tree.m << std::endl;

    // Write true tree to file
    std::ofstream trueTreeFile(OUT_DIR + "/true_tree.nwk");
    if (!trueTreeFile.is_open()) {
        throw std::runtime_error("Failed to open file: " + OUT_DIR + "/true_tree.nwk");
    }
    trueTreeFile << treeToNewick(tree) << std::endl;
    trueTreeFile.close();

    // Write nt to csv file
    std::ofstream trueNtFile(OUT_DIR + "/true_nt.csv");
    if (!trueNtFile.is_open()) {
        throw std::runtime_error("Failed to open file: " + OUT_DIR + "/true_nt.csv");
    }
    for (size_t i = 0; i < tree.nt.size(); i++) {
        trueNtFile << tree.nt[i];
        if (i < tree.nt.size() - 1) {
            trueNtFile << ",";
        }
    }
    trueNtFile << std::endl;
    trueNtFile.close();

    // Initialize inferred tree output
    std::ofstream treeFile(OUT_DIR + "/trees.nexus");
    if (!treeFile.is_open()) {
        throw std::runtime_error("Failed to open file: " + OUT_DIR + "/trees.nexus");
    }
    treeFile << "#NEXUS" << std::endl;
    treeFile << "begin trees;" << std::endl;

    // And nt output
    std::ofstream ntFile(OUT_DIR + "/nt.csv");
    if (!ntFile.is_open()) {
        throw std::runtime_error("Failed to open file: " + OUT_DIR + "/nt.csv");
    }

    // Extract perfect phylogeny
    PerfectPhylo perfectPhylo = extractPerfectPhylo(tree, true); // true means rooted

    State state = constructInitialState(perfectPhylo, NT_INIT);

    //printState(state);
    for(int globalIter = 0; globalIter < N_GLOBAL_ITERS; globalIter++) {

        for(int treeIter = 0; treeIter < N_TREE_ITERS; treeIter++) {
            updateTree(state, rng);
        }

        if(INFER_NT) {
            updateNt(state, rng);
        }

        if(globalIter % SAMPLE_EVERY == 0) {
            Tree treeCurr = stateToTree(state, perfectPhylo, rng);
            // Write tree
            treeFile << "tree STATE_" << globalIter << " = " << treeToNewick(treeCurr) << std::endl;

            // Write nt
            for (size_t i = 0; i < treeCurr.nt.size(); i++) {
                ntFile << treeCurr.nt[i];
                if (i < treeCurr.nt.size() - 1) {
                    ntFile << ",";
                }
            }

            ntFile << std::endl;

        }
        
        // Update progress on same line
        int percentDone = (int)(100.0 * (globalIter + 1) / N_GLOBAL_ITERS);
        std::cout << "\rProgress: " << percentDone << "%" << std::flush;
    }
    std::cout << std::endl;

    // Complete the nexus format
    treeFile << "end;" << std::endl;

    // Close the files
    treeFile.close();
    ntFile.close();

    std::cout << "Success" << std::endl;

    return 0;
}

