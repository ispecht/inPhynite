#pragma once

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include "matrix.h"
#include "tree_inference.h"
#include "newton.h"

// Sample effective population size
void updateNt(
    State& state,
    std::mt19937& rng,
    TridiagonalMatrix& sInv,
    bool forceAccept
); 