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

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [options]" << std::endl;
    std::cout << "Required options:" << std::endl;
    std::cout << "  --nt_init <value>        Initial effective population size" << std::endl;
    std::cout << "  --nt_log_slope <value>   Log slope for nt" << std::endl;
    std::cout << "  --n <value>              Number of leaves" << std::endl;
    std::cout << "  --rooted <0|1>           Root is fixed (0=false, 1=true)" << std::endl;
    std::cout << "  --n_global_iters <value> Number of global iterations" << std::endl;
    std::cout << "  --n_tree_iters <value>   Number of tree iterations" << std::endl;
    std::cout << "  --sample_every <value>   Sample every N iterations" << std::endl;
    std::cout << "  --infer_nt <0|1>         Infer nt (0=false, 1=true)" << std::endl;
    std::cout << "  --out_dir <path>         Output directory" << std::endl;
    std::cout << "Optional:" << std::endl;
    std::cout << "  --seed <value>           Random seed (default: 123)" << std::endl;
}

int main(int argc, char* argv[]) {

    // Parse command line arguments
    std::map<std::string, std::string> args;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg.substr(0, 2) == "--") {
            std::string key = arg.substr(2);
            if (i + 1 < argc) {
                std::string value = argv[i + 1];
                // Check if next arg is a value (not a flag) - allows negative numbers
                if (value.substr(0, 2) != "--") {
                    args[key] = value;
                    i++;
                } else {
                    std::cerr << "Error: Missing value for --" << key << std::endl;
                    printUsage(argv[0]);
                    return 1;
                }
            } else {
                std::cerr << "Error: Missing value for --" << key << std::endl;
                printUsage(argv[0]);
                return 1;
            }
        }
    }

    // Check required arguments
    std::vector<std::string> required = {"nt_init", "nt_log_slope", "n", "rooted", "n_global_iters", 
                                          "n_tree_iters", "sample_every", "infer_nt", "out_dir"};
    for (const auto& req : required) {
        if (args.find(req) == args.end()) {
            std::cerr << "Error: Missing required argument --" << req << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    // Parse arguments
    double NT_INIT = std::stod(args["nt_init"]);
    double NT_LOG_SLOPE = std::stod(args["nt_log_slope"]);
    int N = std::stoi(args["n"]);
    bool ROOTED = (std::stoi(args["rooted"]) != 0);
    int N_GLOBAL_ITERS = std::stoi(args["n_global_iters"]);
    int N_TREE_ITERS = std::stoi(args["n_tree_iters"]);
    int SAMPLE_EVERY = std::stoi(args["sample_every"]);
    bool INFER_NT = (std::stoi(args["infer_nt"]) != 0);
    std::string OUT_DIR = args["out_dir"];
    int SEED = (args.find("seed") != args.end()) ? std::stoi(args["seed"]) : 123;

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
    PerfectPhylo perfectPhylo = extractPerfectPhylo(tree, ROOTED); // true means rooted

    // Write fasta
    std::vector<std::string> sequences = perfectPhyloToSequences(perfectPhylo, 100000); // Length 100000 sequences

    std::ofstream fastaFile(OUT_DIR + "/sequences.fasta");
    if (!fastaFile.is_open()) {
        throw std::runtime_error("Failed to open file: " + OUT_DIR + "/sequences.fasta");
    }
    for (size_t i = 0; i < sequences.size(); i++) {
        fastaFile << ">sequence_" << i << std::endl;
        fastaFile << sequences[i] << std::endl;
    }
    fastaFile.close();


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

