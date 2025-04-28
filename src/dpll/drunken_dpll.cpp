
#include "dpll.h"

drunken_dpll::drunken_dpll(const sat_instance &instance) : dpll(instance) {}

drunken_dpll::~drunken_dpll() {
    delete this->gen;
    delete this->variables;
    delete this->clauses;
}

bool drunken_dpll::unit_elimination(const int &variable, std::vector<int> &current_assignment, std::vector<std::vector<int>> &current_clauses) {
    std::stack<int> unit_variables;
    unit_variables.push(variable);
    while (!unit_variables.empty()) {
        int selected_variable = unit_variables.top();
        unit_variables.pop();

        if (std::find(current_assignment.begin(), current_assignment.end(), -selected_variable) != current_assignment.end()) { 
            return false; 
        }
        if (std::find(current_assignment.begin(), current_assignment.end(), selected_variable) != current_assignment.end()) { continue; }

        current_assignment.push_back(selected_variable);

        for (auto clause_it = current_clauses.begin(); clause_it != current_clauses.end(); ) {
            bool clause_satisfied = false;
            for (auto lit_it = clause_it->begin(); lit_it != clause_it->end(); ) {
                if (*lit_it == selected_variable) {
                    clause_satisfied = true;
                    break;
                }
                else if (*lit_it == -selected_variable) { 
                    lit_it = clause_it->erase(lit_it); 
                }
                else { ++lit_it; }
            }
            if (clause_satisfied) { 
                clause_it = current_clauses.erase(clause_it); 
            } else if (clause_it->empty()) { 
                return false; 
            }
            else if (clause_it->size() == 1) {
                unit_variables.push(clause_it->front()); // new unit clause found
                ++clause_it;
            } else { ++clause_it; }
        }
    }
    return true;
}

int drunken_dpll::random_unassigned_variable(const std::vector<int> &current_assignment) {
    std::vector<int> available_assignments;
    for (const int &variable : *this->variables) {
        bool found = false;
        for (const int &assigned_variable : current_assignment) {
            if (variable == std::abs(assigned_variable)) { 
                found = true; 
                continue; 
            }
        }
        if (found) continue;
        else available_assignments.push_back(variable);
    }
    if (available_assignments.size() == 0) return 0;

    // Make the selection of the variable even more random
    std::uniform_int_distribution<> distrib(0, available_assignments.size() - 1);
    int randomIndex = distrib(*this->gen);
    return available_assignments[randomIndex];
}

bool drunken_dpll::complete(const std::vector<int> &current_assignment) {
    if (current_assignment.size() == this->variables->size()) {
        std::stringstream oss;
        for (const int &variable : current_assignment) {
            oss << std::abs(variable) << (variable > 0 ? " True " : " False ");
        }
        this->assignment = oss.str();
        return true;
    } else return false;
}

bool drunken_dpll::recursive_solver(const int &variable, std::vector<int> current_assignment, std::vector<std::vector<int>> current_clauses) {
    bool unit_elimination_result = this->unit_elimination(variable, current_assignment, current_clauses);
    if (!unit_elimination_result) return false;
    bool complete = this->complete(current_assignment);
    if (complete) return true;
    int unassigned_variable = this->random_unassigned_variable(current_assignment);
    std::uniform_int_distribution<> distrib(0, 1);
    int value = distrib(*this->gen) ? 1 : -1;
    bool flipped_value = distrib(*this->gen);
    complete = this->recursive_solver(flipped_value ? -value * unassigned_variable : value * unassigned_variable, current_assignment, current_clauses);
    if (complete) return true;
    complete = this->recursive_solver(flipped_value ? value * unassigned_variable : -value * unassigned_variable, current_assignment, current_clauses);
    if (complete) return true;
    return false;
}

void drunken_dpll::solve() {
    if (this->debugging) std::cout << "drunken_dpll solver: starting.." << std::endl;

    std::vector<int> assignment;
    this->assignment = "none";

    int unassigned_variable = this->random_unassigned_variable({});
    std::uniform_int_distribution<> distrib(0, 1);
    int value = distrib(*this->gen) ? 1 : -1;

    bool complete = this->recursive_solver(value * unassigned_variable, {}, *this->clauses);
    if (!complete) {
        complete = this->recursive_solver(-value * unassigned_variable, {}, *this->clauses);
    }

    if (this->debugging) {
        std::cout << "\tassignment: " << this->assignment << std::endl;
        std::cout << "\tdrunken_dpll solver: completed" << std::endl;
    }
}

// random restarts

int luby(int i) {
    int k;
    for (k = 1; (1 << (k - 1)) <= i; k++);
    k--;
    if (i == (1 << k) - 1) {
        return (1 << (k - 1));
    }
    return luby(i - (1 << k) + 1);
}

bool drunken_dpll::recursive_solver_with_limit(const int &variable, std::vector<int> current_assignment, std::vector<std::vector<int>> current_clauses, int &decisions, int decision_limit) {
    bool unit_elimination_result = this->unit_elimination(variable, current_assignment, current_clauses);
    if (!unit_elimination_result) return false;

    bool complete = this->complete(current_assignment);
    if (complete) return true;

    if (decisions >= decision_limit) {
        return false;
    }

    decisions++;

    int unassigned_variable = this->random_unassigned_variable(current_assignment);
    std::uniform_int_distribution<> distrib(0, 1);
    int value = distrib(*this->gen) ? 1 : -1;
    bool flipped_value = distrib(*this->gen);

    complete = this->recursive_solver_with_limit(
        flipped_value ? -value * unassigned_variable : value * unassigned_variable,
        current_assignment, current_clauses, decisions, decision_limit
    );
    if (complete) return true;

    complete = this->recursive_solver_with_limit(
        flipped_value ? value * unassigned_variable : -value * unassigned_variable,
        current_assignment, current_clauses, decisions, decision_limit
    );
    if (complete) return true;

    return false;
}

void drunken_dpll::solve_random() {

    if (this->debugging) std::cout << "drunken_dpll solver: starting.." << std::endl;

    int restart_count = 0;
    int decision_limit = luby(1);
    int decisions = 0;

    std::vector<int> assignment;
    this->assignment = "none";

    while (true) {
        // Restart: clear assignment
        assignment.clear();

        std::vector<std::vector<int>> clauses_copy = *this->clauses;

        int unassigned_variable = this->random_unassigned_variable({});
        std::uniform_int_distribution<> distrib(0,1);
        int value = distrib(*this->gen) ? 1 : -1;

        bool complete = this->recursive_solver_with_limit(value * unassigned_variable, assignment, clauses_copy, decisions, decision_limit);

        if (complete) {
            if (this->debugging) {
                std::cout << "\tassignment: " << this->assignment << std::endl;
                std::cout << "\tdrunken_dpll solver: completed" << std::endl;
            }
            return;
        }

        // Restart logic
        restart_count++;
        decision_limit = luby(restart_count + 1);
        if (this->debugging) {
            std::cout << "\tRestarting drunken_dpll after " << decisions << " decisions. New limit: " << decision_limit << std::endl;
        }
        decisions = 0;
    }
}
