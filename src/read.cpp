#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <set>
#include <unordered_set>
#include "sim.h"

std::pair<std::vector<std::string>, std::vector<std::string>> read_fasta(const std::string& filename) {
    std::vector<std::string> names;
    std::vector<std::string> sequences;
    
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    std::string line;
    std::string current_seq = "";
    
    while (std::getline(infile, line)) {
        // Remove trailing whitespace/carriage returns
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' ')) {
            line.pop_back();
        }
        
        if (line.empty()) {
            continue; // Skip empty lines
        }
        
        if (line[0] == '>') {
            // New sequence header
            if (!current_seq.empty()) {
                // Save previous sequence
                sequences.push_back(current_seq);
                current_seq = "";
            }
            // Extract name (skip the '>' character)
            names.push_back(line.substr(1));
        } else {
            // Sequence data - convert to uppercase
            std::transform(line.begin(), line.end(), line.begin(),
                         [](unsigned char c) { return std::toupper(c); });
            current_seq += line;
        }
    }
    
    // Don't forget the last sequence
    if (!current_seq.empty()) {
        sequences.push_back(current_seq);
    }
    
    infile.close();
    
    // Verify we have matching numbers of names and sequences
    if (names.size() != sequences.size()) {
        throw std::runtime_error("Mismatched number of sequence names and sequences");
    }
    
    return {names, sequences};
}

// Filter to positions that have exactly two nucleotides, no Ns or gaps etc
void filter_biallelic_sites(std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return;
    }
    
    size_t seq_length = sequences[0].size();
    
    // Check all sequences have the same length
    for (const auto& seq : sequences) {
        if (seq.size() != seq_length) {
            throw std::runtime_error("All sequences must have the same length");
        }
    }
    
    // Identify which positions to keep
    std::vector<bool> keep_position(seq_length, false);

    int nKeep = 0;
    
    for (size_t pos = 0; pos < seq_length; ++pos) {
        std::unordered_set<char> unique_bases;
        bool has_invalid = false;
        
        // Collect unique bases at this position
        for (const auto& seq : sequences) {
            char base = seq[pos];
            
            // Check if base is valid (A, C, G, or T)
            if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
                has_invalid = true;
                break;
            }
            
            unique_bases.insert(base);
        }
        
        // Keep position only if exactly 2 unique bases and no invalid characters
        if (!has_invalid && unique_bases.size() == 2) {
            keep_position[pos] = true;
            nKeep++;
        }
    }
    
    // Filter each sequence to keep only the valid positions
    for (auto& seq : sequences) {
        std::string filtered_seq;
        filtered_seq.reserve(seq_length); // Reserve space for efficiency
        
        for (size_t pos = 0; pos < seq_length; ++pos) {
            if (keep_position[pos]) {
                filtered_seq += seq[pos];
            }
        }
        
        seq = filtered_seq;
    }

    std::cout << nKeep << " of " << seq_length << " sites are biallelic and have no missing data. Discarding the others." << std::endl;
    std::cout << std::endl;
}


// Filter to positions that have exactly two nucleotides, no Ns or gaps etc
bool remove_incompatible_site(std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return true;
    }
    
    size_t seq_length = sequences[0].size();
    
    // Check all sequences have the same length
    for (const auto& seq : sequences) {
        if (seq.size() != seq_length) {
            throw std::runtime_error("All sequences must have the same length");
        }
    }
    
    // Number of sites each site is incompatible with
    std::vector<int> n_incompatible(seq_length, 0);
    
    for (size_t pos_i = 0; pos_i < seq_length - 1; ++pos_i) {
        for(size_t pos_j = pos_i + 1; pos_j < seq_length; ++pos_j) {
            std::set<std::pair<char, char>> combos;

            for (const auto& seq : sequences) {
                combos.insert({seq[pos_i], seq[pos_j]});
            }

            if(combos.size() > 4) {
                throw std::runtime_error("Should have max 4 combos per site");
            }

            if(combos.size() == 4) {
                // i and j incompatible with one another
                n_incompatible[pos_i]++;
                n_incompatible[pos_j]++;
            }
        }
    }

    // Check if compatible, i.e. all entries are 0
    bool all_zeros = std::all_of(n_incompatible.begin(), n_incompatible.end(), [](int x) { return x == 0; });

    if(all_zeros) {
        return true; // Done!
    }

    auto it = std::max_element(n_incompatible.begin(), n_incompatible.end());
    size_t max_index = std::distance(n_incompatible.begin(), it);
    
    // Remove the most incompatible position
    for (auto& seq : sequences) {
        seq.erase(seq.begin() + max_index);
    }

    return false;
}

void make_compatible(std::vector<std::string>& sequences) {
    size_t seq_length = sequences[0].size();
    int nRemoved = 0;
    std::cout << "Masking sites incompatible with the infinite-sites mutation model using a greedy approach." << std::endl;
    std::cout << "Each site mask costs O(N * L^2), where N is the number of sequences and L is the genome length." << std::endl;
    std::cout << "Please allow sufficient time for this step to run." << std::endl;
    std::cout << std::endl;
    while(!remove_incompatible_site(sequences)){
        nRemoved++;
    }
    std::cout << seq_length - nRemoved << " of " << seq_length << " remaining sites are compatible with the infinite sites model. Discarding the others." << std::endl;
}

std::vector<std::vector<bool>> make_binary(const std::vector<std::string>& sequences) {

    std::vector<std::vector<bool>> bin_seqs(sequences.size(), std::vector<bool>(sequences[0].size(), false));

    for(size_t i = 1; i < sequences.size(); i++) {
        for(size_t pos = 0; pos < sequences[0].size(); pos++) {
            bin_seqs[i][pos] = (sequences[i][pos] != sequences[0][pos]);
        }
    }

    return bin_seqs;
}




// Find the split that is consistent with perfect phylogeny
void buildPerfectPhyloRec(
    bool rooted,
    const std::vector<std::vector<bool>>& sequences,
    const std::vector<bool>& root,
    const std::vector<bool>& isDownstream,
    int nDownstream,
    int idx,
    int& nextAvailIdx,
    std::vector<bool>& splitSoFar,
    PerfectPhylo& soFar
) {

    if(idx >= soFar.w.size()){
        throw std::runtime_error("invalid idx");
    }

    size_t L = sequences[0].size();
    size_t N = sequences.size();

    // Figure out the mutation that appears in the downstream cases the most. This is the one we want to apply!
    std::vector<int> counts(L, 0);
    for(size_t pos = 0; pos < L; pos++) {
        if(splitSoFar[pos]) {
            continue; // already got a mutation at this position
        }

        for(size_t i = 0; i < N; i++) {
            if(isDownstream[i]) {
                if(sequences[i][pos]) {
                    counts[pos]++;
                }
            }
        }
    }

    bool all_zeros = std::all_of(counts.begin(), counts.end(), [](int x) { return x == 0; });
    if(all_zeros) {

        for(size_t i = 0; i < N; i++) {
            if((i==0) && rooted) {
                continue; // Don't count root as leaf
            }

            if(isDownstream[i]) {
                soFar.w[idx]++;
                soFar.leaves[idx].push_back(i);
            }
        }
        return; // Done!
    }

    

    auto it = std::max_element(counts.begin(), counts.end());
    size_t best_pos = std::distance(counts.begin(), it);

    splitSoFar[best_pos] = true;
    

    std::vector<bool> isDownstreamWith = isDownstream; // Has mutation
    std::vector<bool> isDownstreamWithout = isDownstream; // Doesn't have mutation

    int nDownstreamWith = nDownstream;
    int nDownstreamWithout = nDownstream;

    for(size_t i = 0; i < N; i++) {
        if(isDownstream[i]) {
            if(sequences[i][best_pos]) {
                // i has the mutation
                isDownstreamWithout[i] = false;
                nDownstreamWithout--;
            }else{
                // i doesn't have the mutation
                isDownstreamWith[i] = false;
                nDownstreamWith--;
            }
        }
    }

    if(nDownstream != nDownstreamWith + nDownstreamWithout) {
        throw std::runtime_error("splitting incorrect");
    }

    // Process the part of the tree with mutation
    if(nDownstreamWith > 0) {
        std::vector<bool> rootWith = root;

        if(rootWith[best_pos]){
            throw std::runtime_error("should start as false");
        }

        rootWith[best_pos] = true;
        int withIdx = nextAvailIdx;

        soFar.neighbors[withIdx].push_back(idx);
        soFar.neighbors[idx].push_back(withIdx);

        nextAvailIdx++;

        buildPerfectPhyloRec(
            rooted,
            sequences,
            rootWith,
            isDownstreamWith,
            nDownstreamWith,
            withIdx,
            nextAvailIdx,
            splitSoFar,
            soFar
        );

        
    }

    // Process the part of the tree without mutation
    if(nDownstreamWithout > 0) {

        buildPerfectPhyloRec(
            rooted,
            sequences,
            root,
            isDownstreamWithout,
            nDownstreamWithout,
            idx,
            nextAvailIdx,
            splitSoFar,
            soFar
        );
    }
}


PerfectPhylo sequencesToPerfectPhylo(
    const std::vector<std::vector<bool>>& sequences,
    bool rooted
) {

    PerfectPhylo perfectPhylo;
    perfectPhylo.n = sequences.size();
    if(rooted){
        perfectPhylo.n--;
    }
    perfectPhylo.m = sequences[0].size();
    perfectPhylo.w.resize(perfectPhylo.m + 1, 0);
    perfectPhylo.neighbors.resize(perfectPhylo.m + 1, {});
    perfectPhylo.leaves.resize(perfectPhylo.m + 1, {});

    int nextAvailIdx = 1;
    std::vector<bool> splitSoFar(sequences[0].size(), false);

    buildPerfectPhyloRec(
        rooted,
        sequences,
        sequences[0],
        std::vector<bool>(sequences.size(), true),
        sequences.size(),
        0, 
        nextAvailIdx,
        splitSoFar,
        perfectPhylo
    );

    return perfectPhylo;
}




