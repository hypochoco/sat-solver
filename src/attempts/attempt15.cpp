
#include <iostream>
#include <chrono>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <random>
#include <algorithm>
#include <set>
#include <stack>

// cdcl not working that well... if at all...

struct bit {
    int variable = 0;
    bool negative = false;
    int i, j = 0;
    uint64_t bits = 0;

    bit() {}
    bit(const int &variable) {
        this->variable = variable;
        negative = variable < 0;
        i = abs(variable)/64;
        j = abs(variable)%64;
        bits = 1ULL << j;
    }

    int to_variable() const{
        return variable;
    }

    bit neg() const{
        bit neg_bit = *this;
        neg_bit.variable *= -1;
        neg_bit.negative = !negative;
        return neg_bit;
    }

    bool operator==(const bit &other) const {
        return 
        variable == other.variable;
    }
};

struct bitset {
    std::vector<uint64_t> pos_bits;
    std::vector<uint64_t> neg_bits;

    bitset(const int& N) : pos_bits(N/64+1, 0), neg_bits(N/64+1, 0) {}
    bitset(const int& N, const std::unordered_set<int> &nums) : pos_bits(N/64+1, 0), neg_bits(N/64+1, 0) {
        for (const int &num : nums) set(num);
    }
    bitset(const bitset &other) { // copy constructor
        pos_bits = std::vector<uint64_t>(other.pos_bits);
        neg_bits = std::vector<uint64_t>(other.neg_bits);
    }

    bitset emptied() const{
        bitset empty_copy = *this;
        for (uint64_t &bit : empty_copy.pos_bits) { bit = 0; }
        for (uint64_t &bit : empty_copy.neg_bits) { bit = 0; }
        return empty_copy;
    }

    void set(const int &num) {
        if (num >= 0) {
            pos_bits[abs(num)/64] |= 1ULL << (abs(num)%64);
        } else {
            neg_bits[abs(num)/64] |= 1ULL << (abs(num)%64);
        }
    }

    void set (const bit &variable) {
        if (!variable.negative) { // positive
            pos_bits[variable.i] |= variable.bits;
        } else {
            neg_bits[variable.i] |= variable.bits;
        }
    }

    void set (const bitset &other) {
        for (int i=0; i<std::min(size(), other.size()); i++) {
            pos_bits[i] |= other.pos_bits[i];
            neg_bits[i] |= other.neg_bits[i];
        }
    }

    void remove(const bit &variable) {
        if (!variable.negative) { // positive
            if ((pos_bits[variable.i] & variable.bits) > 0) pos_bits[variable.i] ^= variable.bits;
        } else {
            if ((neg_bits[variable.i] & variable.bits) > 0) neg_bits[variable.i] ^= variable.bits;
        }
    }

    void remove_or(const bit &variable) {
        if ((pos_bits[variable.i] & variable.bits) > 0) pos_bits[variable.i] ^= variable.bits;
        if ((neg_bits[variable.i] & variable.bits) > 0) neg_bits[variable.i] ^= variable.bits;
    }

    bool test(const bit &bit) const{
        if (!bit.negative) { // positive
            return (pos_bits[bit.i] & bit.bits) > 0;
        } else { // negative
            return (neg_bits[bit.i] & bit.bits) > 0;
        }
    }

    bool test_or(const bit &bit) const{
        return (pos_bits[bit.i] & bit.bits) > 0 || (neg_bits[bit.i] & bit.bits) > 0;
    }

    bool contains(const bitset &other) const{
        // returns if this bitset contains the other
        if (other.pos_bits.size() > pos_bits.size()) return false; // early stop
        if (other.neg_bits.size() > neg_bits.size()) return false; // early stop

        for (int i=0; i<pos_bits.size(); i++){ // pos bits
            // other.pos_bits cannot have a bit that pos_bits does not
            if ((pos_bits[i] & other.pos_bits[i]) != other.pos_bits[i]) return false;
        }
        for (int i=0; i<neg_bits.size(); i++){ // neg bits
            // other.pos_bits cannot have a bit that pos_bits does not
            if ((neg_bits[i] & other.neg_bits[i]) != other.neg_bits[i]) return false;
        }

        return true;
    }

    int size() const{ return std::max(pos_bits.size(), neg_bits.size()); }

    int num_variables() const{
        int count = 0;
        for (int i=0; i<std::max(pos_bits.size(), neg_bits.size()); i++) {
            if (pos_bits[i] == 0 && neg_bits[i] == 0) continue;
            for (int j=0; j<64; j++) {
                if (pos_bits[i] & (1ULL << j)) count++;
                if (neg_bits[i] & (1ULL << j)) count++;
            }
        }
        return count;
    }

    bitset neg() const{
        bitset negated = *this;
        negated.pos_bits = std::vector<uint64_t>(neg_bits);
        negated.neg_bits = std::vector<uint64_t>(pos_bits);
        return negated;
    }

    bool empty() {
        for (int i=0; i<std::max(pos_bits.size(), neg_bits.size()); i++) {
            if (pos_bits[i] > 0 || neg_bits[i] > 0) return false;
        }
        return true;
    }

    std::vector<int> to_variables() const{
        std::vector<int> variables;
        for (int i=0; i<std::max(pos_bits.size(), neg_bits.size()); i++) {
            if (pos_bits[i] == 0 && neg_bits[i] == 0) continue;
            for (int j=0; j<64; j++) {
                if (pos_bits[i] & (1ULL << j)) variables.push_back(i*64+j);
                if (neg_bits[i] & (1ULL << j)) variables.push_back(-i*64-j);
            }
        }
        return variables;
    }

    std::string to_string() const{
        std::stringstream oss;
        for (int i=0; i<std::max(pos_bits.size(), neg_bits.size()); i++) {
            if (pos_bits[i] == 0 && neg_bits[i] == 0) continue;
            for (int j=0; j<64; j++) {
                if (pos_bits[i] & (1ULL << j)) oss << i*64+j << ' ';
                if (neg_bits[i] & (1ULL << j)) oss << -i*64-j << ' ';
            }
        }
        return oss.str();
    }

    bitset operator&(const bitset &other) const{
        bitset copy = *this;
        for (int i=0; i<other.size(); i++) {
            copy.pos_bits[i] &= other.pos_bits[i];
            copy.neg_bits[i] &= other.neg_bits[i];
        }
        return copy;
    }

    bitset operator&=(const bitset &other) {
        for (int i=0; i<other.size(); i++) {
            pos_bits[i] &= other.pos_bits[i];
            neg_bits[i] &= other.neg_bits[i];
        }
        return *this;
    }
};


struct visitations {
    // absolute basic visitation
    double c = 0;
    std::vector<bitset> *visited;

    visitations() {
        visited = new std::vector<bitset>();
    }

    ~visitations() {
        delete visited;
    }

    void insert(const bitset &assignment) {
        c++;
        for (auto v = visited->begin(); v != visited->end();) {
            if (v->contains(assignment)) {
                v = visited->erase(v);
            } else if (assignment.contains(*v)) {
                return;
            } else ++v;
        }
        visited->push_back(assignment);
    }

    bool test(const bitset &assignment) const{
        for (const bitset &v : *visited) {
            if (assignment.contains(v)) {
                return true;
            }
        }
        return false;
    }

    double count() { return c; }
};

struct sat_problem {
    int N = 0; // number of variables
    std::queue<bit> *unit_queue;
    std::vector<bitset> *clauses;
    std::vector<bitset> *learned_clauses;
    std::unordered_map<int, std::vector<std::pair<int, bit>>> *watch_list;

    sat_problem(const int &_N) : N(_N) {
        // heap allocations
        unit_queue = new std::queue<bit>();
        clauses = new std::vector<bitset>();
        learned_clauses = new std::vector<bitset>();
        watch_list = new std::unordered_map<int, std::vector<std::pair<int, bit>>>();
    }

    ~sat_problem() {
        // delete allocations
        delete unit_queue;
        delete clauses;
        delete learned_clauses;
        delete watch_list;
    }

    void add_clause(const std::unordered_set<int> &clause) {
        if (clause.size() == 0) return; // ignore empty
        else if (clause.size() == 1) { // add to unit queue
            bit num(*clause.begin());
            unit_queue->push(num);
        } else { // add to watchlist 
            bit num1(*clause.begin()); bit num2(*++clause.begin());
            (*watch_list)[num1.to_variable()].push_back({clauses->size(), num2});
            (*watch_list)[num2.to_variable()].push_back({clauses->size(), num1});
            clauses->push_back(bitset(N, clause));
        }
    }

    void add_clause(const bitset &clause) {
        int n = clause.num_variables();
        if (n == 0) return; // ignore empty
        else if (n == 1) { // add to unit queue
            unit_queue->push(clause.to_variables()[0]);
        } else { // add to watchlist 
            std::vector<int> variables = clause.to_variables();
            bit num1(variables[0]); bit num2(variables[1]);
            (*watch_list)[num1.to_variable()].push_back({clauses->size(), num2});
            (*watch_list)[num2.to_variable()].push_back({clauses->size(), num1});
            clauses->push_back(clause);
        }
    }

    void add_learned_clause(const bitset &clause) {
        add_clause(clause);
        learned_clauses->push_back(clause);
    }

    int num_variables() const{ return N; }

    bit pop() {
        if (unit_queue->size() == 0) throw std::runtime_error("cannot pop with queue size 0");
        bit variable = unit_queue->front();
        unit_queue->pop();
        return variable;
    }

    void emplace(const bit &num) { unit_queue->emplace(num); }

    void clear_queue() {
        while (unit_queue->size() > 0) {
            unit_queue->pop();
        }
    }

    int unit_queue_size() { return unit_queue->size(); }

    std::string watch_list_to_string(const std::string &prefix) {
        // variable, other, clause
        std::ostringstream oss;
        for (auto &item : *watch_list) {
            oss << prefix << "variable: " << item.first << std::endl;
            for (auto &pair : item.second) {
                oss << prefix << "\tother: " << pair.second.to_variable() << ", clause: " << (*clauses)[pair.first].to_string() << std::endl;
            }
        }
        return oss.str();
    }
};

bool verify(const bitset &assignment, const sat_problem &sat_data) {
    for (const bitset &clause : *sat_data.clauses) {
        bool satisfied = false;
        for (int i=0; i<assignment.size(); i++) {
            if ((assignment.pos_bits[i] & clause.pos_bits[i]) > 0 || (assignment.neg_bits[i] & clause.neg_bits[i]) > 0) {
                satisfied = true;
                break;
            }
        }
        if (!satisfied) return false;
    }
    return true;
}

int luby(int k) {
    // starts at 1

    int power = 1;
    while (power <= k + 1) power *= 2;
    power /= 2;
    if (power - 1 == k) return power / 2;
    return luby(k - (power - 1));
}

struct vsids {

    // variables
    int N; // num variables
    int c; int j; float c_add; float c_decay;
    std::unordered_map<int, float>* variable_scores; // variable -> value
    std::set<std::pair<float, int>>* sorted_variable_scores; // value -> variable

    // normalization
    float normalization_value;
    float max_score_threshold;

    // perturbation
    std::mt19937 gen;

    vsids(const int &_N) : N(_N) {
        // note
            // neither of these two make a large difference
            // having vsids makes a decent difference though

        // variables
        c = 0;
        // j = 10; c_add = 1; c_decay = 0.95;
        j = 1000; c_add = 1; c_decay = 0.5;
        variable_scores = new std::unordered_map<int, float>();
        sorted_variable_scores = new std::set<std::pair<float, int>>();

        // normalization
        normalization_value = 1e2;
        max_score_threshold = 1e3;

        // perturbations
        gen = std::mt19937(std::random_device{}());

        // vairable initialization
        for (int i=0; i<N; i++) {
            (*variable_scores)[i+1] = 0;
            sorted_variable_scores->insert({0, i+1});
        }
    }

    ~vsids() {
        delete variable_scores;
        delete sorted_variable_scores;
    }

    void normalize() {
        // normalize map and set
        std::set<std::pair<float, int>> sorted_variable_scores_copy;
        for (int i=0; i< N; i++) {
            (*variable_scores)[i+1] /= normalization_value;
            sorted_variable_scores_copy.insert({(*variable_scores)[i+1], i+1});
        }
        *sorted_variable_scores = std::move(sorted_variable_scores_copy);
    }

    void add(const bitset &clause) {
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

    void decay(const bitset &clause) {
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

    int next(const bitset &assignment) const{
        // find the next variable that hasn't been assigned
        for (auto it = sorted_variable_scores->rbegin(); it != sorted_variable_scores->rend(); ++it) {
            bit variable(it->second);
            if (!assignment.test_or(variable)) { // not assigned
                return variable.to_variable();
            }
        }
        return -1; // no more variables to be assigned
    }

    void perturb() {
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

};

struct cdcl {

    // data structures
        // decision level (int) -> variable (int), clause (index, -1 if a decision)
    std::vector<std::pair<bit,int>> *trail; 

    cdcl() { 
        
        trail = new std::vector<std::pair<bit,int>>(); 
    }

    ~cdcl() { delete trail; }

    void implication(const bit &variable, const int &clause_index) {
        // add a decision to the graphs
        trail->push_back({variable, clause_index});
    }

    void decision(const bit &variable) {
        // add an implication to the graph
        implication(variable, -1);
    }

    std::pair<bitset, bitset> resolve_conflict(const bitset &conflict_clause, const std::vector<bitset> *clauses) {
        // iterate through the graph to construct a learned clause

        // cache variables at the current decision level
        bitset variables_at_decision_level = conflict_clause.emptied();
        for (auto it = trail->rbegin(); it != trail->rend(); ++it) {
            variables_at_decision_level.set(it->first);
            if (it->second == -1) break;
        }

        // create learned clause
        auto it = trail->rbegin();
        bitset learned_clause(conflict_clause);
        while ((learned_clause.neg() & variables_at_decision_level).num_variables() > 1) {

            // trace backward through the trail
            if (it->second == -1) {
                if ((learned_clause.neg() & variables_at_decision_level).num_variables() > 1) {
                    throw std::runtime_error("failed!!");
                } else break;
            }
            if (!learned_clause.test(it->first.neg())) { // skip if not in current learned clause
                ++it; continue;
            } 
            bitset reason_clause = (*clauses)[it->second]; // reason clause

            // update learned clause
            learned_clause.set(reason_clause);
            learned_clause.remove_or(it->first);

            // next iteration
            ++it;
        }

        // find second highest decision level
        int inverse_decision_level = 0;
        for (auto it = trail->rbegin(); it != trail->rend(); ++it) {
            if (it->second == -1) inverse_decision_level++;
            if (learned_clause.test_or(it->first)) {
                if (it->second != -1) inverse_decision_level++;
                break;
            }
        }

        // remove variables greater than the second highest decision level
        int i = 0;
        bitset backtrack_variables = conflict_clause.emptied();
        for (auto it = trail->rbegin(); it != trail->rend(); ++it) {
            backtrack_variables.set(it->first);
            if (it->second == -1) i++;
            if (i == inverse_decision_level) {
                ++it;
                trail->resize(trail->size() - std::distance(trail->rbegin(), it));
                break;
            }
        }

        // return learned clause and backtrack variables
        return {learned_clause, backtrack_variables};
    }

    std::string trail_to_string() {
        std::stringstream oss;
        for (auto it = trail->begin(); it != trail->end(); ++it) {
            if (it->second == -1) oss << it->first.to_variable() << "d ";
            else oss << it->first.to_variable() << "i ";
        }
        return oss.str();
    }

};

struct sat_solver {
    // random number generation
    std::mt19937 gen;

    // solver data structures
    sat_problem &sat_data;
    visitations &visited;

    // luby restarts
    int base;
    int iteration_count;
    int max_iterations;

    // vsids
    vsids decision;

    // cdcl
    cdcl backtrack;

    sat_solver(sat_problem &_sat_data, visitations &_visited, const int &N) : sat_data(_sat_data), visited(_visited), decision(N) {

        // random generation
        gen = std::mt19937(std::random_device{}());

        // luby restarts, default values
        base = 1; iteration_count = 0; max_iterations = 1;
    }

    // dpll functions

    std::optional<bitset> unit_propagation(bitset &assignment) {
        // return conflict 

        // handle learned clauses
        bitset assignment_overlap = assignment.emptied();
        for (const bitset &learned_clause : *sat_data.learned_clauses) {
            assignment_overlap.set(assignment & learned_clause.neg());
        }
        for (const int &variable : assignment_overlap.to_variables()) {
            sat_data.emplace(bit(variable));
        }

        // iterate over variables marked for deletion
        while (sat_data.unit_queue_size() > 0) {

            // pop item from queue
            const bit &bit_num = sat_data.pop();
            const bit &bit_num_neg = bit_num.neg();

            // assign variable in assignment, cache negative
            assignment.set(bit_num);

            // update watched literals
            std::vector<std::tuple<int,int,bit>> add_list;
            for (auto pair = (*sat_data.watch_list)[bit_num_neg.to_variable()].begin(); pair != (*sat_data.watch_list)[bit_num_neg.to_variable()].end();) {
                const int bitset_clause_index = pair->first;
                const bitset bitset_clause = (*sat_data.clauses)[bitset_clause_index];
                const bit bit_other = pair->second;

                // iterate over variables in clause
                bool new_watched_literal_found = false;
                for (const int &num : bitset_clause.to_variables()) { // note: could do randomization here
                    const bit bit_new_num(num);

                    // check bit conditions (not other, negative not assigned)
                    if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num.neg())) {
                        new_watched_literal_found = true;
                        add_list.push_back({bit_new_num.to_variable(),bitset_clause_index,bit_other});

                        for (std::pair<int, bit> &other_pair : (*sat_data.watch_list)[bit_other.to_variable()]) { // update other
                            if (other_pair.first == bitset_clause_index && other_pair.second == bit_num_neg) {
                                other_pair.second = bit_new_num;
                                break;
                            }
                        }
                        break;
                    }
                }

                // no new watched literal found
                if (!new_watched_literal_found) {
                    if (assignment.test(bit_other.neg())) {
                        // vsids update
                        decision.add(bitset_clause);
                        decision.decay(bitset_clause);

                        // additions
                        for (const std::tuple<int,int,bit> &addition : add_list) {
                            (*sat_data.watch_list)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
                        }
                        return bitset_clause;
                    }
                    else if (!assignment.test_or(bit_other)) {
                        if (assignment.test(bit_other.neg())) { return bitset_clause; } // other neg already assigned
                        backtrack.implication(bit_other, bitset_clause_index);
                        sat_data.emplace(bit_other);
                    }
                    ++pair;
                } else { // watched literal found
                    pair = (*sat_data.watch_list)[bit_num_neg.to_variable()].erase(pair);
                }
            }
            // additions
            for (const std::tuple<int,int,bit> &addition : add_list) {
                (*sat_data.watch_list)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
            }
        }
        return std::nullopt;
    }

    std::string solve(bitset &assignment) {

        // variables
        std::stack<std::pair<bitset,bit>> stack; // allocate this onto the heap??
        std::uniform_int_distribution<> dis(0, 1); // random polarity

        // initial iteration
        const std::optional<bitset> init_conflict = unit_propagation(assignment);
        if (init_conflict.has_value()) { return "unsat"; }
        int unassigned_variable = decision.next(assignment); // vsids variable selection
        if (unassigned_variable == -1) return "sat"; // stopping condition
        int polarity = dis(gen) * 2 - 1;
        stack.push({assignment, bit(polarity*unassigned_variable)});
        stack.push({assignment, bit(-polarity*unassigned_variable)});

        // cdcl cache original
        bitset original_assignment(assignment);

        // solving
        while (stack.size() > 0) {

            // note
                // without random restarts and without luby without cdcl
                // C1597_024.cnf, 859.59, SAT passed

            // // luby restarts
            // if (iteration_count > max_iterations) return "stopped";

            // pop assignment from stack, prepare iteration
            std::pair<bitset,bit> pair = stack.top(); stack.pop();
            assignment = pair.first; bit unassigned_bit = pair.second;
            sat_data.clear_queue();
            sat_data.emplace(unassigned_bit);

            // cdcl decision
            backtrack.decision(unassigned_bit);

            // std::cout << "visited count: " << visited.count() << std::endl;

            // // visited
            // bitset visited_assignment(assignment);
            // visited_assignment.set(unassigned_bit);
            // if (visited.test(visited_assignment)) continue;

            // // debugging
            // std::cout << std::endl;
            // std::cout << "assignment before unit prop: " << assignment.to_string() << ", variable: " << unassigned_bit.to_variable() << std::endl;

            // unit propagation
            const std::optional<bitset> conflict = unit_propagation(assignment);
            if (conflict) { // conflict found

                // // visited assignments
                // visited.insert(visited_assignment);






                // // we need visitations... 
                // // there's no way of telling sometimes if we've already visited something... 

                // // code feels very spaghetti like now... 
                // // what do we do about that... oofs...





                // // // luby restarts
                // // iteration_count++;

                // // // debugging
                // // std::cout << "conflict hit" << std::endl;

                // // cdcl resolve conflict (learned clause, backtrack variables)
                // std::pair<bitset,bitset> backtrack_result = backtrack.resolve_conflict(conflict.value(), sat_data.clauses);
                // bitset learned_clause = backtrack_result.first;
                // bitset backtrack_variables = backtrack_result.second;

                // // // debugging
                // // std::cout << "current assignment: " << assignment.to_string() << std::endl;
                // // std::cout << "learned clause: " << learned_clause.to_string() << std::endl;
                // // std::cout << "backtrack variables: " << backtrack_variables.to_string() << std::endl;

                // // std::cout << "stack: " << std::endl;
                // // std::stack<std::pair<bitset,bit>> stack_copy(stack);
                // // while (!stack_copy.empty()) {
                // //     auto pair = stack_copy.top();
                // //     stack_copy.pop();
                // //     std::cout << "\tbit: " << pair.second.to_variable() << ", assignment: " << pair.first.to_string() << std::endl; 
                // // }

                // // stop unit clauses
                // if (learned_clause.num_variables() == 1) throw std::runtime_error("stop learned unit clause");

                // sat_data.add_learned_clause(learned_clause);

                // // remove from stack
                // while (stack.size() > 0) {
                //     std::pair<bitset,bit> pair = stack.top();
                //     if (!(pair.first & backtrack_variables).empty() || backtrack_variables.test_or(pair.second)) {
                //         stack.pop();
                //     } else break;
                // }      
                
                // // std::cout << "cut stack: " << std::endl;
                // // std::stack<std::pair<bitset,bit>> stack_cut_copy(stack);
                // // while (!stack_cut_copy.empty()) {
                // //     auto pair = stack_cut_copy.top();
                // //     stack_cut_copy.pop();
                // //     std::cout << "\tbit: " << pair.second.to_variable() << ", assignment: " << pair.first.to_string() << std::endl; 
                // // }

                // // proposed new addition
                // // std::cout << "proposed new addition: " << std::endl;
                // if (!stack.empty()) {
                //     std::pair<bitset, bit> other = stack.top();
                //     stack.push({other.first, other.second.neg()});
                //     // std::cout << "\tbit: " << other.second.neg().to_variable() << ", assignment: " << other.first.to_string() << std::endl; 
                // }
                
                continue;
                // throw std::runtime_error("stop");

            }

            // assign variable
            int unassigned_variable = decision.next(assignment); // vsids variable selection
            if (unassigned_variable == -1) return "sat"; // stopping condition
            int polarity = dis(gen) * 2 - 1;
            stack.push({assignment, bit(polarity*unassigned_variable)});
            stack.push({assignment, bit(-polarity*unassigned_variable)});
        }
        return "unsat";
    }

    std::string restart(bitset &assignment) {
        
        // notes:
            // visitations and restarts (parallelization) -> something doesn't feel right
            // cdcl

            // check memory usage?
        
        // // luby restarts
        // base = 1e5;
        // for (int i=1; i<1e6; i++) {
        //     iteration_count = 0;
        //     max_iterations = luby(i) * base;
        //     bitset restart_assignment(assignment);
        //     decision.perturb();
        //     std::string result = solve(restart_assignment);
        //     if (result == "sat") { // found result
        //         assignment = restart_assignment;
        //         return result;
        //     } else if (result == "unsat") return result; // found result
        // }
        // return "stopped";

        // no luby restarts
        return solve(assignment);
    }
};

struct unit_tests {

    static bool test_bitset_contains() {

        bitset A(150);
        A.set(1);
        A.set(50);
        A.set(125);
        A.set(-125);

        bitset B(150);
        B.set(1);
        B.set(50);
        B.set(-125);

        if (A.contains(B) != true) return false;
        if (B.contains(A) != false) return false;

        return true;
    }

    static bool test_bitset_and() {

        bitset A(150);
        A.set(1); A.set(50); A.set(-125);
        A.set(125);

        bitset B(150);
        B.set(1); B.set(50); B.set(-125);
        B.set(149);

        bitset C(150);
        C.set(1); C.set(50); C.set(-125);

        if (C.contains((A&B)) && (A&B).contains(C)) return true;
        else return false; 
    }

    static bool test_bitset_emptied() {

        bitset A(150);
        A.set(1);
        A.set(50);
        A.set(125);
        A.set(-125);

        if (A.emptied().num_variables() == 0) return true;
        else return false;
    }

    static void run_tests() {
        std::string test_bit_contains_status = test_bitset_contains()? "passed" : "failed";
        std::cout << "test_bitset_contains: " << test_bit_contains_status << std::endl;

        std::string test_bit_and_status = test_bitset_and()? "passed" : "failed";
        std::cout << "test_bitset_and: " << test_bit_and_status << std::endl;

        std::string test_bit_emptied_status = test_bitset_emptied()? "passed" : "failed";
        std::cout << "test_bitset_emptied: " << test_bit_emptied_status << std::endl;
    }
};

int main(int argc, char* argv[]) {

    // // debugging -- unit tests
    // unit_tests::run_tests();
    // return 0;
    
    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    // parse input

    // input file error handling
    if (argc != 2) { std::cerr << "Usage: " << argv[0] << " <input.cnf>" << std::endl; return 1; }
    std::ifstream cnfFile(argv[1]);
    if (!cnfFile) { std::cerr << "Error: Cannot open file " << argv[1] << std::endl; return 1; }

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

    // watch list
        // variable -> list (clause, other variable)
    
    sat_problem sat_data(N);
    
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
        sat_data.add_clause(unorderded_set_clause); // add clause to data structure
    }

    // implementation

    // debugging -- logging start
    std::cout << "starting solve" << std::endl;

    bitset assignment(N);
    visitations visisted;
    sat_solver sat_solver(sat_data, visisted, N);
    std::string result = sat_solver.restart(assignment);
    
    // log results
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    std::cout << "total time: " << time.count() << std::endl;
    std::cout << "result: " << result << std::endl;
    if (result == "sat") {
        std::cout << "assignment: " << assignment.to_string() << std::endl;
        std::string verification = verify(assignment, sat_data)? "verified" : "failed";
        std::cout << "verification: " << verification << std::endl;
    }

    // // final output
    // auto end = std::chrono::high_resolution_clock::now();
    // const std::chrono::duration<double> time = end - start;
    // std::string formatted_result = "";
    // if (result == "sat") formatted_result = "SAT";
    // else if (result == "unsat") formatted_result = "UNSAT";
    // else if (result == "stopped") formatted_result = "STOPPED";
    // size_t pos = std::string(argv[1]).find_last_of("/\\");
    // std::string filename = (pos != std::string::npos) ? std::string(argv[1]).substr(pos + 1) : "";
    // std::cout << "{\"Instance\": \"" << filename << "\", \"Time\": " << std::setprecision(2) << std::fixed << time.count() << ", \"Result\": \"" << formatted_result << "\"";
    // if (result == "sat") {
    //     std::string verification = verify(assignment, sat_data)? "Passed" : "Failed";
    //     std::cout << ", \"Verification\": \"" << verification << "\"";
    // }
    // std::cout << "}" << std::endl;

    // default condition
    return 0;
}
