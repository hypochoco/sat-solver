
#include <iostream>
#include <chrono>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <random>
#include <algorithm>

struct bit {
    int variable = 0;
    bool negative = false;
    int i, j = 0;
    uint64_t bits = 0;

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

    void set(const int &num) {
        if (num >= 0) {
            pos_bits[abs(num)/64] |= 1ULL << (abs(num)%64);
        } else {
            neg_bits[abs(num)/64] |= 1ULL << (abs(num)%64);
        }
    }

    void set (const bit &variable ) {
        if (!variable.negative) { // positive
            pos_bits[variable.i] |= variable.bits;
        } else {
            neg_bits[variable.i] |= variable.bits;
        }
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

    int size() const{ return std::max(pos_bits.size(), neg_bits.size()); }

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
            bit num1(*clause.begin());
            bit num2(*++clause.begin());
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

bool partial_verify(const bitset &assignment, const sat_problem &sat_data, const bool &verbose = false) {

    // ignores clauses that do not have all variables assigned
    // otherwise, performs regular verification

    std::vector<int> assigned_variables = assignment.to_variables();
    if (assigned_variables.size() == 0) return true; // early stop, no variables assigned

    // debugging -- assignemnt variables
    if (verbose) std::cout << "[partial verify] assignment variables: " << assignment.to_string() << std::endl;

    for (const bitset &clause : *sat_data.clauses) {
        // check all variables in clause are in the assignment
        std::vector<int> clause_variables = clause.to_variables();
        bool all_variables_assigned = true;
        for (const int &clause_variable : clause_variables) {
            bit bit_clause_variable(clause_variable);
            if (!assignment.test_or(bit_clause_variable)) {
                all_variables_assigned = false;

                // // debugging -- early stop
                // if (verbose) std::cout << "\tcould not find variable " << clause_variable << " in assignment" << std::endl;

                break;
            }
        }
        if (!all_variables_assigned) continue;

        // // debugging -- checking variables
        // if (verbose) std::cout << "[partial verify] all variables found, attempting verification" << std::endl;

        // check correct assignment
        bool satisfied = false;
        for (int i=0; i<assignment.size(); i++) {
            if ((assignment.pos_bits[i] & clause.pos_bits[i]) > 0 || (assignment.neg_bits[i] & clause.neg_bits[i]) > 0) {
                satisfied = true;
                break;
            }
        }
        if (!satisfied) {

            // debugging -- failing conditions
            if (verbose) std::cout << "could not satisfy clause: " << clause.to_string() << std::endl;

            return false;
        }
    }
    return true;
}

struct sat_solver {
    // random number generation
    std::mt19937 gen;

    // solver data structures
    sat_problem &sat_data;

    sat_solver(sat_problem &_sat_data) : sat_data(_sat_data) {
        gen = std::mt19937(std::random_device{}());
    }

    // dpll functions

    bool unit_propagation(bitset &assignment) {
        // return conflict 

        // iterate over variables marked for deletion
        while (sat_data.unit_queue_size() > 0) {

            // pop item from queue
            const bit &bit_num = sat_data.pop();
            const bit &bit_num_neg = bit_num.neg();
            if (assignment.test(bit_num_neg)) return false; // early stop

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
                for (const int &num : bitset_clause.to_variables()) { 
                    const bit bit_new_num(num);

                    // check bit conditions (not other, negative not assigned)
                    if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num.neg())) {
                        new_watched_literal_found = true;
                        add_list.push_back({bit_new_num.to_variable(),bitset_clause_index,bit_other});

                        // debugging -- add list
                        std::cout << "adding to add list: " << bitset_clause_index << ", new add list size: " << add_list.size() << std::endl;

                        for (std::pair<int, bit> &other_pair : (*sat_data.watch_list)[bit_other.to_variable()]) { // update other
                            if (other_pair.first == bitset_clause_index && other_pair.second == bit_num_neg) {
                                other_pair.second = bit_new_num;
                                break;
                            }
                        }
                        break;
                    }
                }

                std::cout << "add list size 1: " << add_list.size() << std::endl;

                // no new watched literal found
                if (!new_watched_literal_found) {
                    if (assignment.test(bit_other.neg())) {

                        // additions
                        for (const std::tuple<int,int,bit> &addition : add_list) {

                            // debugging -- adding
                            std::cout << "adding clause: " << std::get<1>(addition) << ": " << (*sat_data.clauses)[std::get<1>(addition)].to_string() << std::endl;

                            (*sat_data.watch_list)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
                        }

                        return false;
                    }
                    else if (!assignment.test_or(bit_other)) sat_data.emplace(bit_other);
                    ++pair;
                } else { // watched literal found

                    // debugging -- erase
                    std::cout << "removing clause: " << pair->first << ": " << (*sat_data.clauses)[pair->first].to_string() << std::endl;

                    pair = (*sat_data.watch_list)[bit_num_neg.to_variable()].erase(pair);
                }

                std::cout << "add list size 2: " << add_list.size() << std::endl;


            }
            // additions
            for (const std::tuple<int,int,bit> &addition : add_list) {

                // debugging -- adding
                std::cout << "adding clause: " << std::get<1>(addition) << ": " << (*sat_data.clauses)[std::get<1>(addition)].to_string() << std::endl;

                (*sat_data.watch_list)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
            }

            std::cout << std::endl;
        }
        return true;
    }

    std::string solve(bitset &assignment) {
        // recursive solver

        // visitations TODO
        
        // unit propagation
        if (!unit_propagation(assignment)) return "unsat";

        // degbugging -- partial verification
        bool partial_verification = partial_verify(assignment, sat_data);
        if (!partial_verification) {
            std::cout << "partial verification failed" << std::endl;
            partial_verify(assignment, sat_data, true);

            // debugging the watchlist
            std::cout << "watchlist:" << std::endl;
            std::cout << sat_data.watch_list_to_string("\t") << std::endl;

            throw std::runtime_error("partial verification failed");
        }

        // randomized variable selection
        std::vector<int> unassigned_variables;
        for (int i=1; i<=sat_data.num_variables(); i++) {
            if (!assignment.test_or(i)) { unassigned_variables.push_back(i); }
        }
        int unassigned_variable = -1;
        if (unassigned_variables.size() > 0) {
            std::uniform_int_distribution<> dis(0, unassigned_variables.size() - 1);
            unassigned_variable = unassigned_variables[dis(gen)];
        }

        // stopping condition
        if (unassigned_variable == -1) return "sat";

        // branching
        bitset assignment_true = assignment;
        bit unassigned_bit_true(abs(unassigned_variable));
        sat_data.emplace(unassigned_bit_true);
        std::string result_true = solve(assignment_true);
        if (result_true == "sat") {
            assignment = assignment_true;
            return "sat";
        }
        sat_data.clear_queue(); // clear queue
        bitset assignment_false = assignment;
        bit unassigned_bit_false(-abs(unassigned_variable));
        sat_data.emplace(unassigned_bit_false);
        std::string result_false = solve(assignment_false);
        if (result_false == "sat") {
            assignment = assignment_false;
            return "sat";
        }
        
        return "unsat";
    }

    std::string restart(bitset &assignment) {
        return solve(assignment);
    }
};

struct unit_tests {
    static void run_tests() {}
};

int main(int argc, char* argv[]) {

    // debugging -- unit tests
    unit_tests::run_tests();
    
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
    // visitations visited;
    
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

            // unorderded_set_clause.insert(int_value);
        }
        sat_data.add_clause(unorderded_set_clause); // add clause to data structure
    }

    // debugging -- init data structures
    std::cout << "init watchlist:" << std::endl;
    std::cout << sat_data.watch_list_to_string("\t") << std::endl;

        // note:
            // it looks like things are disappearing from the watch list
            // where does it go... 

    // implementation
    bitset assignment(N);
    sat_solver sat_solver(sat_data);
    std::string result = sat_solver.restart(assignment);
    
    // debugging -- final time
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    std::cout << "total time: " << time.count() << std::endl;

    // log results
    std::cout << "result: " << result << std::endl;
    if (result == "sat") {
        std::cout << "assignment: " << assignment.to_string() << std::endl;

        // partial verification
        std::string partial_verification = partial_verify(assignment, sat_data)? "verified" : "failed";
        std::cout << "partial verification: " << partial_verification << std::endl;

        // normal verification
        std::string verification = verify(assignment, sat_data)? "verified" : "failed";
        std::cout << "verification: " << verification << std::endl;
    }

    // default condition
    return 0;

    // notes
        // returning sat when it should be unsat
        // failing verification on things that are sat
            // C1597_024

        // toy test 2 returns sat, but it should be unsat
            // verification fails as well
            // partial verifications to find when this breaks...
}
