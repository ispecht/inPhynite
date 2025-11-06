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
#include "read.h"

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [options]" << std::endl;
    std::cout << "Required options:" << std::endl;
    std::cout << "  --rooted <0|1>           Root is fixed (0=false, 1=true)" << std::endl;
    std::cout << "  --n_global_iters <value> Number of global iterations" << std::endl;
    std::cout << "  --n_tree_iters <value>   Number of tree iterations" << std::endl;
    std::cout << "  --sample_every <value>   Sample every N iterations" << std::endl;
    std::cout << "  --infer_nt <0|1>         Infer nt (0=false, 1=true)" << std::endl;
    std::cout << "  --out_dir <path>         Output directory" << std::endl;
    std::cout << "Optional:" << std::endl;
    std::cout << "  --seed <value>           Random seed (default: 123)" << std::endl;
    std::cout << "  --nt_init <value>        Initial effective population size" << std::endl;
    std::cout << "  --nt_log_slope <value>   Log slope for nt" << std::endl;
    std::cout << "  --n <value>              Number of leaves" << std::endl;
    std::cout << "  --in_file <path>         Input FASTA file" << std::endl;
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
    std::vector<std::string> required = {"rooted", "n_global_iters", 
                                          "n_tree_iters", "sample_every", "infer_nt", "out_dir"};
    for (const auto& req : required) {
        if (args.find(req) == args.end()) {
            std::cerr << "Error: Missing required argument --" << req << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    // Parse arguments
    bool ROOTED = (std::stoi(args["rooted"]) != 0);
    int N_GLOBAL_ITERS = std::stoi(args["n_global_iters"]);
    int N_TREE_ITERS = std::stoi(args["n_tree_iters"]);
    int SAMPLE_EVERY = std::stoi(args["sample_every"]);
    bool INFER_NT = (std::stoi(args["infer_nt"]) != 0);
    std::string OUT_DIR = args["out_dir"];

    // Optional args
    int SEED = (args.find("seed") != args.end()) ? std::stoi(args["seed"]) : 123;
    double NT_INIT = (args.find("nt_init") != args.end()) ? std::stod(args["nt_init"]) : -1001.0; // Default
    double NT_LOG_SLOPE = (args.find("nt_log_slope") != args.end()) ? std::stod(args["nt_log_slope"]) : -1001.0;
    int N = (args.find("n") != args.end()) ? std::stoi(args["n"]) : -1;

    std::string IN_FILE = (args.find("in_file") != args.end()) ? args["in_file"] : "";

    PerfectPhylo perfectPhylo; // Initialize data in the form of perfect phylo

    std::mt19937 rng;
    rng.seed(SEED);

    // If there's an input file, read it
    if(IN_FILE != "") {

        NT_INIT = 1.0; // Starting value of effective population size inference

        std::pair<std::vector<std::string>, std::vector<std::string>> fasta = read_fasta(IN_FILE);
        std::vector<std::string> names = fasta.first;
        std::vector<std::string> in_sequences = fasta.second;
        filter_biallelic_sites(in_sequences);
        make_compatible(in_sequences);
        std::vector<std::vector<bool>> bin_seqs = make_binary(in_sequences);
        perfectPhylo = sequencesToPerfectPhylo(bin_seqs, false);

    }else{

        if(N < 0) {
            throw std::runtime_error("For simulated data, n (number of leaves) must be specified");
        }

        if(NT_INIT < -1000.0) {
            throw std::runtime_error("For simulated data, nt_init (initial effective population size) must be specified");
        }

        if(NT_LOG_SLOPE < -1000.0) {
            throw std::runtime_error("For simulated data, nt_log_slope (slope of log effective population size) must be specified");
        }

        Tree tree = simulateTree(
            N,
            NT_INIT,
            NT_LOG_SLOPE,
            rng
        );

        std::cout << "Number of mutations on tree: " << tree.m << std::endl;

        // Extract perfect phylogeny
        perfectPhylo = extractPerfectPhylo(tree, ROOTED); // true means rooted


        // Write true tree to file
        std::ofstream trueTreeFile(OUT_DIR + "/true_tree.nwk");
        if (!trueTreeFile.is_open()) {
            throw std::runtime_error("Failed to open file: " + OUT_DIR + "/true_tree.nwk");
        }
        trueTreeFile << treeToNewick(tree) << std::endl;
        trueTreeFile.close();

        // Write true nt to csv file
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
    }

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


    State state = constructInitialState(perfectPhylo, NT_INIT);

    TridiagonalMatrix sInv(
        std::vector<double>(state.n-1, 0.0),
        std::vector<double>(state.n-2, 0.0)
    );

    if(INFER_NT) {
        updateNt(state, rng, sInv, true); // Force accept 1st change in case of terrible initial guess
    }

    //printState(state);

    // Time the Gibbs loop
    auto start = std::chrono::high_resolution_clock::now();

    for(int globalIter = 0; globalIter < N_GLOBAL_ITERS; globalIter++) {

        for(int treeIter = 0; treeIter < N_TREE_ITERS; treeIter++) {
            updateTree(
                state,
                treeIter % (state.m + state.n - 2),
                rng
            );
        }

        if(INFER_NT) {
            updateNt(state, rng, sInv, false);
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
        //std::cout << state.nt[0] << std::endl;
    }

    std::cout << std::endl;


    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Sampling took: " << (double) duration.count() / (1000000.0) << " seconds" << std::endl;




    std::cout << std::endl;

    // Complete the nexus format
    treeFile << "end;" << std::endl;

    // Close the files
    treeFile.close();
    ntFile.close();

    std::cout << "Success" << std::endl;

    return 0;
}

