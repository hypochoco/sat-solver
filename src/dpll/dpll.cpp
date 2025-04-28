
#include "dpll.h"

dpll::dpll(const sat_instance &instance) {

    // debugging
    if (this->debugging) std::cout << "dpll constructor: building internal representation.." << std::endl;

    // random
    std::random_device rd;
    this->gen = new std::mt19937(rd());

    // internal representation
    this->variables = new std::unordered_set<int>();
    this->clauses = new std::vector<std::vector<int>>();
    for (const auto &line : *instance.clauses) {
        std::vector<int> clause; // note: does not account for duplicates
        std::string value;
        std::stringstream ss(line);
        while (ss >> value) {
            int int_value = std::stoi(value);
            if (int_value == 0) break;
            clause.push_back(int_value);
            variables->emplace(std::abs(int_value));
        }
        this->clauses->push_back(clause);
    }

    // debugging
    if (this->debugging) {
        std::ostringstream oss;
        for (const int &variable : *variables) {
            oss << variable << " ";
        }
        std::cout << "\tvariables: " << oss.str() << std::endl;
        std::cout << "\tdpll constructor: completed" << std::endl;
    }
}

dpll::~dpll() {
    delete this->gen;
    delete this->variables;
    delete this->clauses;
}

bool dpll::unit_elimination(const int &variable, std::vector<int> &current_assignment, std::vector<std::vector<int>> &current_clauses) {
    std::stack<int> unit_variables;
    unit_variables.push(variable);
    while (!unit_variables.empty()) {
        int selected_variable = unit_variables.top();
        unit_variables.pop();
        if (std::find(current_assignment.begin(), current_assignment.end(), -selected_variable) != current_assignment.end()) { return false; }
        if (std::find(current_assignment.begin(), current_assignment.end(), selected_variable) != current_assignment.end()) { continue; }
        current_assignment.push_back(selected_variable);
        for (auto clause_it = current_clauses.begin(); clause_it != current_clauses.end(); ) {
            bool clause_satisfied = false;
            for (auto lit_it = clause_it->begin(); lit_it != clause_it->end(); ) {
                if (*lit_it == selected_variable) {
                    clause_satisfied = true;
                    break;
                }
                else if (*lit_it == -selected_variable) { lit_it = clause_it->erase(lit_it); }
                else { ++lit_it; }
            }
            if (clause_satisfied) { clause_it = current_clauses.erase(clause_it); }
            else if (clause_it->empty()) { return false; }
            else if (clause_it->size() == 1) {
                unit_variables.push(clause_it->front()); // new unit clause found
                ++clause_it;
            } else { ++clause_it; }
        }
    }
    return true;
}

int dpll::random_unassigned_variable(const std::vector<int> &current_assignment) {
    // return random unassigned variable
    std::vector<int> available_assignments;
    for (const int &variable : *this->variables) {
        bool found = false;
        for (const int &assigned_variable : current_assignment) {
            if (variable == std::abs(assigned_variable)) { found = true; continue; }
        }
        if (found) continue;
        else available_assignments.push_back(variable);
    }
    if (available_assignments.size() == 0) return 0;
    std::uniform_int_distribution<> distrib(0, available_assignments.size() - 1);
    int randomIndex = distrib(*this->gen);
    return available_assignments[randomIndex];
}

bool dpll::complete(const std::vector<int> &current_assignment) {
    // checks for assignment completion, if complete write to assignmet
    if (current_assignment.size() == this->variables->size()) { // complete, build assignment string
        std::stringstream oss;
        for (const int &variable : current_assignment) {
            oss << std::abs(variable) << (variable>0? " True " : " False ");
        }
        this->assignment = oss.str();
        return true;
    } else return false;
}

bool dpll::recursive_solver(const int &variable, std::vector<int> current_assignment, std::vector<std::vector<int>> current_clauses) {

    // unit propogation
    bool unit_elimination_result = this->unit_elimination(variable, current_assignment, current_clauses);
    if (!unit_elimination_result) return false;

    // check completion
    bool complete = this->complete(current_assignment);
    if (complete) return true;

    // choose random unassigned variable and random value
    int unassigned_variable = this->random_unassigned_variable(current_assignment);
    std::uniform_int_distribution<> distrib(0,1);
    int value = distrib(*this->gen) ? 1 : -1;

    // branch value 0
    complete = this->recursive_solver(value*unassigned_variable, current_assignment, current_clauses);
    if (complete) return true;

    // branch value 1
    complete = this->recursive_solver(-value*unassigned_variable, current_assignment, current_clauses);
    if (complete) return true;

    // no possible assignment
    return false;
}

void dpll::solve() {

    // algo (recursive)
        // eliminate unit literals
        // choose a variable to branch on
        // backtrack when necessary

    // debugging
    if (this->debugging) std::cout << "dpll solver: starting.." << std::endl;

    // data structures and variables
    std::vector<int> assignment;
    this->assignment = "none";

    // solver
    int unassigned_variable = this->random_unassigned_variable({});
    std::uniform_int_distribution<> distrib(0,1);
    int value = distrib(*this->gen) ? 1 : -1;
    bool complete = this->recursive_solver(value*unassigned_variable, {}, *this->clauses);
    if (!complete) {
        complete = this->recursive_solver(-value*unassigned_variable, {}, *this->clauses);
    }

    // debugging 
    if (this->debugging) {
        std::cout << "\tassignment: " << this->assignment << std::endl;
        std::cout << "\tdpll solver: completed" << std::endl;
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

bool dpll::recursive_solver_with_limit(const int &variable, std::vector<int> current_assignment, std::vector<std::vector<int>> current_clauses, int &decisions, int decision_limit) {

    // unit propagation
    bool unit_elimination_result = this->unit_elimination(variable, current_assignment, current_clauses);
    if (!unit_elimination_result) return false;

    // check completion
    bool complete = this->complete(current_assignment);
    if (complete) return true;

    // Check decision limit for restart
    if (decisions >= decision_limit) {
        return false; // Force a restart
    }

    decisions++;

    // choose random unassigned variable and random value
    int unassigned_variable = this->random_unassigned_variable(current_assignment);
    std::uniform_int_distribution<> distrib(0,1);
    int value = distrib(*this->gen) ? 1 : -1;

    // branch value 0
    complete = this->recursive_solver_with_limit(value * unassigned_variable, current_assignment, current_clauses, decisions, decision_limit);
    if (complete) return true;

    // branch value 1
    complete = this->recursive_solver_with_limit(-value * unassigned_variable, current_assignment, current_clauses, decisions, decision_limit);
    if (complete) return true;

    // no possible assignment
    return false;
}

void dpll::solve_random() {

    if (this->debugging) std::cout << "dpll solver: starting.." << std::endl;

    int restart_count = 0;
    int decision_limit = luby(1); // Start with luby(1)
    int decisions = 0;

    std::vector<int> assignment;
    this->assignment = "none";

    while (true) {
        // Start fresh
        assignment.clear();

        std::vector<std::vector<int>> clauses_copy = *this->clauses;

        int unassigned_variable = this->random_unassigned_variable({});
        std::uniform_int_distribution<> distrib(0,1);
        int value = distrib(*this->gen) ? 1 : -1;

        bool complete = this->recursive_solver_with_limit(value * unassigned_variable, assignment, clauses_copy, decisions, decision_limit);

        if (complete) {
            if (this->debugging) {
                std::cout << "\tassignment: " << this->assignment << std::endl;
                std::cout << "\tdpll solver: completed" << std::endl;
            }
            return;
        }

        // Restart logic
        restart_count++;
        decision_limit = luby(restart_count + 1);
        if (this->debugging) {
            std::cout << "\tRestarting solver after " << decisions << " decisions. New limit: " << decision_limit << std::endl;
        }
        decisions = 0;
    }
}
