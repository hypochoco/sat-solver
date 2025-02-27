
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
#include <thread>
#include <future>
#include <limits>

// log16

// data structures

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

    bitset() {}
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

struct sat_problem {
    int N = 0; // number of variables
    std::queue<bit> *unit_queue;
    std::vector<bitset> *clauses;
    std::unordered_map<int, std::vector<std::pair<int, bit>>> *watch_list;

    sat_problem(const int &_N) : N(_N) {
        // heap allocations
        unit_queue = new std::queue<bit>();
        clauses = new std::vector<bitset>();
        watch_list = new std::unordered_map<int, std::vector<std::pair<int, bit>>>();
    }

    sat_problem(const sat_problem& other) : N(other.N) {
        unit_queue = new std::queue<bit>(*other.unit_queue);
        clauses = new std::vector<bitset>(*other.clauses);
        watch_list = new std::unordered_map<int, std::vector<std::pair<int, bit>>>(*other.watch_list);
    }

    ~sat_problem() {
        // delete allocations
        delete unit_queue;
        delete clauses;
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

    ~vsids() {
        delete variable_scores;
        delete sorted_variable_scores;
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

    void normalize() {

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

};

// struct moms {

//     int next(const sat_problem &sat_data, bitset &assignment) const{
//         // find all unsatisfied clauses and min size
//         int min_clause_size = 0;
//         std::vector<int> unsat_clause_indices;
//         for (int i=0; i<(*sat_data.clauses).size(); i++) {
//             if (((*sat_data.clauses)[i] & assignment).empty()) {
//                 unsat_clause_indices.push_back(i);
//                 min_clause_size = std::min(min_clause_size, (*sat_data.clauses)[i].num_variables());
//             }
//         }

//         // completely satisfied
//         if (unsat_clause_indices.size() == 0) return -1;

//         // count variable freq in smallest clauses
//         std::unordered_map<int, int> variable_freq;
//         for (int i : unsat_clause_indices) {
//             if ((*sat_data.clauses)[i].num_variables() == min_clause_size) {
//                 for (const int &variable : (*sat_data.clauses)[i].to_variables()) {
//                     variable_freq[variable]++;
//                 }
//             }
//         }

//         // sort
//         std::set<std::pair<int, int>> sorted_variable_freq;
//         for (int i=0; i<sat_data.num_variables(); i++) {
//             sorted_variable_freq.insert({variable_freq[i+1], i});
//         }

//         // return next variable
//         return sorted_variable_freq.rbegin()->second;
//     }

// };

struct sat_solver {
    // random number generation
    std::mt19937 gen;

    // solver data structures
    sat_problem &sat_data;
    std::stack<std::pair<bitset,bit>> *stack;

    // vsids
    vsids decision;
    // moms decision_moms;

    sat_solver(sat_problem &_sat_data) : sat_data(_sat_data), decision(_sat_data.num_variables()) {
        gen = std::mt19937(std::random_device{}());
        stack = new std::stack<std::pair<bitset,bit>>();
    }

    ~sat_solver() {
        delete stack;
    }

    std::optional<bitset> unit_propagation(bitset &assignment) {
        // return conflict 

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

    std::string solve(bitset &assignment, std::atomic<bool> &value_found) {

        // variables
        std::uniform_int_distribution<> dis(0, 1);

        // initial iteration
        decision.perturb(); // randomize initial decision
        const std::optional<bitset> init_conflict = unit_propagation(assignment);
        if (init_conflict.has_value()) { return "unsat"; }
        int unassigned_variable = decision.next(assignment); // vsids variable selection
        // int unassigned_variable = decision_moms.next(sat_data, assignment); // moms vairable selection
        if (unassigned_variable == -1) return "sat"; // stopping condition
        int polarity = dis(gen) * 2 - 1;
        stack->push({assignment, bit(polarity*unassigned_variable)});
        stack->push({assignment, bit(-polarity*unassigned_variable)});

        // solving
        while (stack->size() > 0) {

            // multithreading stop
            if (value_found.load()) return "unsat";

            // pop assignment from stack, prepare iteration
            std::pair<bitset,bit> pair = stack->top(); stack->pop();
            assignment = pair.first; bit unassigned_bit = pair.second;
            sat_data.clear_queue();
            sat_data.unit_queue->push(unassigned_bit);
            sat_data.emplace(unassigned_bit);

            // unit propagation
            const std::optional<bitset> conflict = unit_propagation(assignment);
            if (conflict) { continue; }

            // assign variable
            int unassigned_variable = decision.next(assignment); // vsids variable selection
            // int unassigned_variable = decision_moms.next(sat_data, assignment); // moms variable selection
            if (unassigned_variable == -1) return "sat"; // stopping condition
            int polarity = dis(gen) * 2 - 1;
            stack->push({assignment, bit(polarity*unassigned_variable)});
            stack->push({assignment, bit(-polarity*unassigned_variable)});
        }
        return "unsat";
    }
};

struct pool {
    std::mutex result_mutex;
    sat_problem &sat_data;
    bitset default_assignment;
    bitset found_assignment;
    std::string result = "unsat";
    std::atomic<bool> value_found{false};
    
    pool(sat_problem &_sat_data) : sat_data(_sat_data) {
        default_assignment = bitset(_sat_data.num_variables());
        found_assignment = bitset(_sat_data.num_variables());
    }

    void task(const std::vector<int> branch) {

        // run solver on branch
        sat_problem sat_data_copy(sat_data);
        for (const int &variable : branch) sat_data_copy.emplace(bit(variable));
        sat_solver solver(sat_data_copy);
        bitset task_assignment(default_assignment);
        std::string task_result = solver.solve(task_assignment, value_found);

        // get and store 
        if (task_result == "sat") {
            std::lock_guard<std::mutex> lock(result_mutex);
            if (!value_found.load()) {
                value_found.store(true);
                result = "sat";
                found_assignment = task_assignment;
            }
        }
    }

    std::vector<std::vector<int>> create_branches(const int &n) {

        // note:
            // dlcs -> not good
            // moms -> it's okay... 

        // min number of variables
        int num_branches = 1 << n;

        // find smallest clause
        int min_clause_size = INT_MAX;
        for (const bitset &clause : *sat_data.clauses) {
            int m = clause.num_variables();
            if (m <= 1 || m < num_branches) continue;
            min_clause_size = std::min(min_clause_size, m);
        }

        // count variable freq in smallest clauses
        std::unordered_map<int, int> variable_freq;
        for (const bitset &clause : *sat_data.clauses) {
            if (clause.num_variables() == min_clause_size) {
                for (const int &variable : clause.to_variables()) {
                    variable_freq[variable]++;
                }
            }
        }

        // sort
        std::set<std::pair<int, int>> sorted_variable_freq;
        for (int i=0; i<sat_data.num_variables(); i++) {
            sorted_variable_freq.insert({variable_freq[i+1], i});
        }

        // create branches
        std::vector<std::vector<int>> branches;
        for (int i = 0; i < num_branches; i++) {
            std::vector<int> b;
            auto it = sorted_variable_freq.rbegin();
            for (int j = 0; j < n; j++) {
                if (it == sorted_variable_freq.rend()) break;
                int sign = (i & (1 << j)) ? -1 : 1;
                b.push_back(sign * it->second);
                ++it;
            }
            branches.push_back(b);
        }

        // debugging -- branches
        std::cout << "branches: " << std::endl;
        for (const std::vector<int> &x : branches) {
            std::cout << "\t";
            for (const int &y : x) std::cout << y << " ";
            std::cout << std::endl;
        }

        // look at the branches and see if there's something off...
        // what about a smarter multithreading strategy...
        // looking at the stack rather than whatever else... 


        // delay branches??
        // put things 


        // different, cooler branching strategy...
        // maybe hold off on branches that take a long while... 
        // when do branches end???


        // print out when those branches end... 


        return branches;
    }

    std::pair<std::string, bitset> run() {

        // note:
            // try out different heuristics...
            // try out different vsids parameters...
                // perturb it every once in a while??
            // try more branches...

            // what is life and stuff in this place...
            // how can we start getting the stuff we can't... 

        // branches
        std::vector<std::vector<int>> branches = create_branches(3);

        // threads
        std::vector<std::thread> threads;
        for (const std::vector<int> b : branches) {
            threads.emplace_back(&pool::task, this, b);
        }
        for (auto& t : threads) { if (t.joinable()) t.join(); }

        // return found output
        return {result, found_assignment};
    }
};

// helper functions

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

sat_problem preprocess(const int &argc, char* argv[]) {
    // parse input

    // input file error handling
    // if (argc != 2) { std::cerr << "Usage: " << argv[0] << " <input.cnf>" << std::endl; return 1; }
    std::ifstream cnfFile(argv[1]);
    // if (!cnfFile) { std::cerr << "Error: Cannot open file " << argv[1] << std::endl; return 1; }

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

    return sat_data;
}

int main(int argc, char* argv[]) {
    
    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    // process inputs
    sat_problem sat_data = preprocess(argc, argv);

    // implementation -- threading
    pool p(sat_data);
    std::pair<std::string, bitset> result_pair = p.run();
    std::string result = result_pair.first;
    bitset assignment = result_pair.second;
    
    // // log results
    // auto end = std::chrono::high_resolution_clock::now();
    // const std::chrono::duration<double> time = end - start;
    // std::cout << "total time: " << time.count() << std::endl;
    // std::cout << "result: " << result << std::endl;
    // if (result == "sat") {
    //     std::cout << "assignment: " << assignment.to_string() << std::endl;
    //     std::string verification = verify(assignment, sat_data)? "verified" : "failed";
    //     std::cout << "verification: " << verification << std::endl;
    // }

    // final output
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    std::string formatted_result = "";
    if (result == "sat") formatted_result = "SAT";
    else if (result == "unsat") formatted_result = "UNSAT";
    else if (result == "stopped") formatted_result = "STOPPED";
    size_t pos = std::string(argv[1]).find_last_of("/\\");
    std::string filename = (pos != std::string::npos) ? std::string(argv[1]).substr(pos + 1) : "";
    std::cout << "{\"Instance\": \"" << filename << "\", \"Time\": " << std::setprecision(2) << std::fixed << time.count() << ", \"Result\": \"" << formatted_result << "\"";
    if (result == "sat") {
        std::string verification = verify(assignment, sat_data)? "Passed" : "Failed";
        std::cout << ", \"Verification\": \"" << verification << "\"";
    }
    std::cout << "}" << std::endl;

    // default condition
    return 0;
}
