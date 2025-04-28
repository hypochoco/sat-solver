
#include <set>
#include <random>
#include <unordered_map>

#include "decision.h"
#include "data_structures.h"


void vsids::perturb() {

    // random number generation
    std::uniform_real_distribution<float> dist(0.0, 1.0);

    // perturb all values in map and set
    std::set<std::pair<float, int>> sorted_variable_scores_copy;
    for (int i=0; i< N; i++) {
        (*variable_scores)[i+1] /= dist(gen);
        sorted_variable_scores_copy.insert({(*variable_scores)[i+1], i+1});
    }
    *sorted_variable_scores = std::move(sorted_variable_scores_copy);
}

void vsids::normalize() {

    // normalize map and set
    std::set<std::pair<float, int>> sorted_variable_scores_copy;
    for (int i=0; i< N; i++) {
        (*variable_scores)[i+1] /= normalization_value;
        sorted_variable_scores_copy.insert({(*variable_scores)[i+1], i+1});
    }
    *sorted_variable_scores = std::move(sorted_variable_scores_copy);

    // perturb after normalization
    perturb();
}

void vsids::add(const bitset &clause) {

    // require normalization
    bool require_normalize = false;

    // give a conflict clause, update all the variables by c_add
    for (const int &variable : clause.to_variables()) {
        int abs_variable = abs(variable);

        // remove old
        sorted_variable_scores->erase({(*variable_scores)[abs_variable], abs_variable});

        // insert new value
        (*variable_scores)[abs_variable] += c_add;
        sorted_variable_scores->insert({(*variable_scores)[abs_variable], abs_variable});

        // check normalization
        if ((*variable_scores)[abs_variable] >= max_score_threshold) require_normalize = true;
    }

    // normalize if threshold passed
    if (require_normalize) normalize();
}

void vsids::decay(const bitset &clause) {

    c++;
    if (c%j == 0) { // update on j iterations
        c = 0;

        // decay clauses
        for (const int &variable : clause.to_variables()) {
            int abs_variable = abs(variable);

            // remove old
            sorted_variable_scores->erase({(*variable_scores)[abs_variable], abs_variable});

            // insert new value
            (*variable_scores)[abs_variable] *= c_decay;
            sorted_variable_scores->insert({(*variable_scores)[abs_variable], abs_variable});
        }
    }
}

int vsids::next(const bitset &assignment) const{
    // find the next variable that hasn't been assigned
    for (auto it = sorted_variable_scores->rbegin(); it != sorted_variable_scores->rend(); ++it) {
        bit variable(it->second);
        if (!assignment.test_or(variable)) { // not assigned
            return variable.to_variable();
        }
    }
    return -1; // no more variables to be assigned
}
