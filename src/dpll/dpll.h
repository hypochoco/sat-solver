
#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <stack>
#include <algorithm>
#include <unordered_set>
#include <random>

#include "../solver.h"

struct dpll : solver {
    dpll(const sat_instance &instance);
    ~dpll();

    virtual bool unit_elimination(const int &variable, std::vector<int> &current_assignment, std::vector<std::vector<int>> &current_clauses);
    virtual int random_unassigned_variable(const std::vector<int> &current_assignment);
    virtual bool complete(const std::vector<int> &current_assignment);
    virtual bool recursive_solver(const int &variable, std::vector<int> current_assignment, std::vector<std::vector<int>> current_clauses);
    virtual void solve();

    virtual bool recursive_solver_with_limit(const int &variable, std::vector<int> current_assignment, std::vector<std::vector<int>> current_clauses, int &decisions, int decision_limit);
    virtual void solve_random();

    bool debugging = false;
    std::mt19937 *gen;
    std::unordered_set<int> *variables;
    std::vector<std::vector<int>> *clauses;
};

struct drunken_dpll : dpll {
    drunken_dpll(const sat_instance &instance);
    ~drunken_dpll();

    virtual bool unit_elimination(const int &variable, std::vector<int> &current_assignment, std::vector<std::vector<int>> &current_clauses);
    virtual int random_unassigned_variable(const std::vector<int> &current_assignment);
    virtual bool complete(const std::vector<int> &current_assignment);
    virtual bool recursive_solver(const int &variable, std::vector<int> current_assignment, std::vector<std::vector<int>> current_clauses);
    virtual void solve();

    virtual bool recursive_solver_with_limit(const int &variable, std::vector<int> current_assignment, std::vector<std::vector<int>> current_clauses, int &decisions, int decision_limit);
    virtual void solve_random();
};
