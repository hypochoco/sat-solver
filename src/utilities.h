
#pragma once

#include <vector>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "data_structures.h"
#include "watch_list.h"


inline bool verify(const bitset &assignment, const std::vector<bitset> &clauses) {
    for (const bitset &clause : clauses) {
        bool satisfied = false;
        for (int i=0; i<assignment.size(); i++) {
            if ((assignment.pos_bits[i] & clause.pos_bits[i]) > 0 || (assignment.neg_bits[i] & clause.neg_bits[i]) > 0) {
                satisfied = true;
                break;
            }
        }
    }
    return true;
}

inline int luby(int k) {
    // starts at 1

    int power = 1;
    while (power <= k + 1) power *= 2;
    power /= 2;
    if (power - 1 == k) return power / 2;
    return luby(k - (power - 1));
}

inline std::tuple<int, watch_list*, std::vector<bitset>*> preprocess(const int &argc, char* argv[]) {
    // number of variables, watch_list, clauses

    // parse input
    std::ifstream cnfFile(argv[1]);

    // initial stats
    int stats[2];
    std::string line;
    while (std::getline(cnfFile, line)) {
        if (line[0] == 'c') continue; // ignore comments
        else if (line[0] == 'p') { // num variables and num clauses
            int count = 0;
            std::string value;
            std::stringstream ss(line);
            while (std::getline(ss, value, ' ')) { // iterate over line
                if (value == "p") continue;
                if (value == "cnf") continue;
                stats[count] = std::stoi(value);
                count++;
            }
            break;
        } else break;
    }
    int N = stats[0]; // number of variables
    int M = stats[1]; // number of clauses

    // data structures
    watch_list *wc = new watch_list(N);
    std::vector<bitset> *clauses = new std::vector<bitset>();
    
    // reduce gaps in variables, create a mapping of variables
    int max_abs_variable = 0;
    int reduction_count = 0;
    std::unordered_map<int, int> reduction_map;

    // iterate over input cnf file
    while (std::getline(cnfFile, line)) {
        if (line[0] == 'c') continue; // ignore comments
        if (line[0] == 'p') continue; // ignore stats

        // initialization
            // unit clauses  -> queued for unit propagation
            // unordered set -> bitmask
            // watch list    -> 2 literals per clause
        std::string value;
        std::stringstream ss(line);
        std::unordered_set<int> unorderded_set_clause;
        while (std::getline(ss, value, ' ')) { // iterate over variable
            if (value == "0") continue;
            int int_value = std::stoi(value);
            max_abs_variable = std::max(max_abs_variable, abs(int_value)); // for logging purposes

            if (reduction_map[abs(int_value)] == 0) { // if default value, new item
                reduction_count++;
                reduction_map[abs(int_value)] = reduction_count;
            }
            if (int_value >= 0) unorderded_set_clause.insert(reduction_map[abs(int_value)]);
            else unorderded_set_clause.insert(-reduction_map[abs(int_value)]);
        }
        // add to watch list and clauses
        if (unorderded_set_clause.size() == 0) continue;
        else if (unorderded_set_clause.size() == 1) {
            wc->add_unit_clause(bit(*unorderded_set_clause.begin()));
        } else {
            wc->add_clause(
                bit(*unorderded_set_clause.begin()), 
                bit(*++unorderded_set_clause.begin()),
                clauses->size()
            );
            clauses->push_back(bitset(N, unorderded_set_clause));
        }
    }

    return {N, wc, clauses};
}