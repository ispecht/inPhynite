#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <string>
#include "helpers.h"

void removeValue(std::vector<int>& vec, int value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    
    if (it == vec.end()) {
        throw std::invalid_argument("Value not found in vector");
    }
    
    vec.erase(it);
}