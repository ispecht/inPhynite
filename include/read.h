#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <map>
#include <random>
#include "sim.h"

std::pair<std::vector<std::string>, std::vector<std::string>> read_fasta(const std::string& filename);
void filter_biallelic_sites(std::vector<std::string>& sequences);
void make_compatible(std::vector<std::string>& sequences);
std::vector<std::vector<bool>> make_binary(const std::vector<std::string>& sequences);

PerfectPhylo sequencesToPerfectPhylo(
    const std::vector<std::vector<bool>>& sequences,
    bool rooted
);