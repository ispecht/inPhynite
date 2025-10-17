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


int main() {

    double ntInit = 1.0;
    double ntLogSlope = 0.0;
    int n = 10;


    const int SEED = 0;
    std::mt19937 rng;
    rng.seed(SEED);

    Tree tree = simulateTree(
        n,
        ntInit,
        ntLogSlope,
        rng
    );

    // printTree(tree);
    PerfectPhylo perfectPhylo = extractPerfectPhylo(tree, true); // true means rooted
    //printPerfectPhylo(perfectPhylo);

    State state = constructInitialState(perfectPhylo, ntInit);

    //printState(state);

    for(int i = 0; i < 1000; i++) {
        updateTree(state, rng);
    }

    std::cout << "Success" << std::endl;

    return 0;
}

