
#include <set>
#include <unordered_map>
#include <iostream>

#include "decision.h"
#include "data_structures.h"


int moms::next(const bitset &assignment) const{

    // note
        // not considering the size after assignments

    // find smallest assignments
    int min_size = INT_MAX;
    std::vector<int> unsat_clauses;
    for (int i=0; i<clauses.size(); i++) {
        if ((clauses[i] & assignment).empty()) {
            int n = clauses[i].num_variables();
            if (n < min_size) {
                min_size = n;
                unsat_clauses = std::vector<int>();
                unsat_clauses.push_back(i);
            } else if (n == min_size) {
                unsat_clauses.push_back(i);
            }
        }
    }

    // check termination
    if (unsat_clauses.size() == 0) {
        if (assignment.num_variables() < N) { // satisfied, but not all vars assigned
            for (int i=1; i<N+1; i++) {
                if (!assignment.test_or(i)) return i;
            }
        } else return -1;
    }

    // find max number of occurrences
    std::unordered_map<int,int> freq_map;
    for (const int &i : unsat_clauses) {
        for (const int &variable : clauses[i].to_variables()) {
            freq_map[abs(variable)]++;
        }
    }
    std::set<std::pair<int,int>> sorted_freq;
    for (const std::pair<int,int> &it : freq_map) sorted_freq.insert({it.second, it.first});

    // ensure not assigned
    for (auto it = sorted_freq.rbegin(); it != sorted_freq.rend(); ++it) {
        bit variable(it->second);
        if (assignment.test_or(variable)) continue;
        return it->second;
    }

    // error
    throw std::runtime_error("moms error!");
}
