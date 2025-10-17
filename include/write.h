#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <random>
#include "sim.h"
#include "tree_inference.h"

Tree stateToTree(
    const State& state,
    PerfectPhylo perfectPhylo, // passing by value because we will update this
    std::mt19937& rng
);

std::string treeToNewick(const Tree& tree);

std::vector<std::string> perfectPhyloToSequences(const PerfectPhylo& perfectPhylo, int L);
