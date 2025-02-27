
#pragma once

#include <queue>
#include <vector>
#include <iostream>

#include "decision.h"
#include "data_structures.h"


struct watch_list {

    int N = 0;
    std::queue<bit> *unit_queue;
    std::unordered_map<int, std::vector<std::pair<int, bit>>> *watched_literals;

    inline watch_list(const int &_N) : N(_N) {
        unit_queue = new std::queue<bit>();
        watched_literals = new std::unordered_map<int, std::vector<std::pair<int, bit>>>();
    }

    inline watch_list(const watch_list &other) : N(other.N) {
        unit_queue = new std::queue<bit>(*other.unit_queue);
        watched_literals = new std::unordered_map<int, std::vector<std::pair<int, bit>>>(*other.watched_literals);
    }

    inline ~watch_list() {
        delete unit_queue;
        delete watched_literals;
    }

    inline void add_clause(const bit &num1, const bit &num2, const int &clause_index) {
        (*watched_literals)[num1.to_variable()].push_back({clause_index, num2});
        (*watched_literals)[num2.to_variable()].push_back({clause_index, num1});
    }
    
    inline void add_unit_clause(const bit &variable) {
        unit_queue->push(variable);
    }

    bool unit_propagation(bitset &assignment, const std::vector<bitset> &clauses, decision &d);
    void calibrate(const bitset &assignment, const std::vector<bitset> &clauses);
    int num_variables() const;
    bit pop();
    void emplace(const bit &num);
    void clear_queue();
    int unit_queue_size() const;
    std::string watch_list_to_string(const std::string &prefix, const std::vector<bitset> *clauses) const;

    void sanity_check(const bitset &assignment, const std::string &debug) const;

};
