#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include "sim.h"


int main() {
    std::cout << "hi" << std::endl;

    const int SEED = 0;
    std::mt19937 rng;
    rng.seed(SEED);

    Tree tree = simulateTree(
        10,
        10.0,
        0.0,
        rng
    );

    //printTree(tree);

    PerfectPhylo perfectPhylo = extractPerfectPhylo(tree, true);
    //printPerfectPhylo(perfectPhylo);

    return 0;
}

