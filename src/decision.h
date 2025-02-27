
#pragma once

#include <set>
#include <random>
#include <unordered_map>

#include "data_structures.h"

struct decision {

    // variables
    int N; // num variables

    inline decision(const int &_N) : N(_N) {}
    inline virtual ~decision() {};

    // vsids functions
    virtual void perturb() = 0;
    virtual void normalize() = 0;
    virtual void add(const bitset &clause) = 0;
    virtual void decay(const bitset &clause) = 0;

    // moms functions

    // decision functions
    virtual int next(const bitset &assignment) const = 0;

};

struct vsids : decision {
    // variables
    int c; // decay count
    int j; // decay iterations
    float c_add; // score increase amount
    float c_decay; // decay multiplier
    std::unordered_map<int, float>* variable_scores; // variable -> value
    std::set<std::pair<float, int>>* sorted_variable_scores; // value -> variable
    float normalization_value; // normalization
    float max_score_threshold; // normalization
    std::mt19937 gen; // randomization

    inline vsids(const int &_N) : decision(_N) {
        // variables
        c = 0;
        // j = 10; c_add = 1; c_decay = 0.95;
        j = 1000; c_add = 1; c_decay = 0.5; // generally faster
        variable_scores = new std::unordered_map<int, float>();
        sorted_variable_scores = new std::set<std::pair<float, int>>();

        // normalization
        normalization_value = 1e2;
        max_score_threshold = 1e3;

        // perturbations
        gen = std::mt19937(std::random_device{}());

        // vairable initialization
        for (int i=0; i<N; i++) {
            (*variable_scores)[i+1] = 1;
            sorted_variable_scores->insert({0, i+1});
        }
    }

    inline ~vsids() {
        delete variable_scores;
        delete sorted_variable_scores;
    }

    void perturb();
    void normalize();
    void add(const bitset &clause);
    void decay(const bitset &clause);
    int next(const bitset &assignment) const;
};

struct moms : decision {
    // variables
    const std::vector<bitset> &clauses;

    inline moms(const int &_N, const std::vector<bitset> &_clauses) : decision(_N), clauses(_clauses) {}
    inline ~moms() {}

    int next(const bitset &assignment) const;

    // unused
    inline void perturb() {};
    inline void normalize() {};
    inline void add(const bitset &clause) {};
    inline void decay(const bitset &clause) {};
};
