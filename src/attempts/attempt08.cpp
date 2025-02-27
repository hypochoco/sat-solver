
#include <iostream>
#include <chrono>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <random>

// notes:
    // better memory management techniques
    // better visualizations of information

// data structures

struct bit {
    int _i = 0; // location in a bits vector (0->N)
    int _j = 0; // location of bit (2*k, k -> [0,32])
    int _variable = 0; // variable
    uint64_t bits = 0; // bits representation

    bit(const int &variable) { set(variable); }

    void set(const int &variable) { // sets the variable
        _variable = variable;
        int value = 2*abs(variable);
        _i = value/64;
        _j = value%64;
        bits |= (1ULL << _j);
        if (variable < 0) bits |= (1ULL << (_j-1));
    }

    bit neg() const{ // returns a negative copy
        bit num_neg = *this;
        num_neg._variable *= -1;
        num_neg.bits ^= 1ULL<<(_j-1);
        return num_neg;
    }

    int to_variable() const { // returns variable representation
        return _variable;
    }

    bool operator==(const bit &other) const { // == operator definition
        return _variable == other._variable;
    }
};

struct bitset {
    int c; // number of variable stored
    std::vector<uint64_t> pos_bits;
    std::vector<uint64_t> neg_bits;

    bitset(const int &_N) : pos_bits((2*_N + 63) / 64, 0), neg_bits((2*_N + 63) / 64, 0) {}
    bitset(const int &_N, const std::unordered_set<int> &clause) : pos_bits((2*_N + 63) / 64, 0), neg_bits((2*_N + 63) / 64, 0) {
        for (const int &variable : clause) set(variable);
    }

    void set(const int &variable) {
        c++;
        int value = 2*abs(variable);
        if (variable >= 0) {
            pos_bits[value/64] |= (1ULL << (value%64));
        } else {
            neg_bits[value/64] |= (1ULL << (value%64));
            neg_bits[(value-1)/64] |= (1ULL << ((value-1)%64));
        }
    }

    void set(const bit &variable) {
        c++;
        if (variable.to_variable() >= 0) {
            pos_bits[variable._i] |= variable.bits;
        } else {
            neg_bits[variable._i] |= variable.bits;
        }
    }

    bool test(const int &variable) const{
        int value = 2*abs(variable);
        if (variable >= 0) {
            if ((pos_bits[value/64] & (1ULL << (value%64))) == (1ULL << (value%64))) return true;
        } else {
            if ((neg_bits[value/64] & (1ULL << (value%64))) == (1ULL << (value%64))) return true;
        }
        return false;
    }

    bool test_or(const int &variable) const{
        int value = 2*abs(variable);
        if ((pos_bits[value/64] & (1ULL << (value%64))) == (1ULL << (value%64))) return true;
        if ((neg_bits[value/64] & (1ULL << (value%64))) == (1ULL << (value%64))) return true;
        return false;
    }

    bool test(const bit &variable) const{
        if (variable.to_variable() >= 0) {
            if ((pos_bits[variable._i] & variable.bits) == variable.bits) return true;
        } else {
            if ((neg_bits[variable._i] & variable.bits) == variable.bits) return true;
        }
        return false;
    }

    bool test_or(const bit &variable) const{
        if ((pos_bits[variable._i] & variable.bits) == variable.bits) return true;
        if ((neg_bits[variable._i] & variable.bits) == variable.bits) return true;
        return false;
    }

    int count() const{ return c; }

    // int size() const{ return N; }
    int size() const{ return pos_bits.size(); }

    std::vector<int> to_variables() const{
        std::vector<int> variables;
        for (int i=0; i<pos_bits.size(); i++) {
            for (int j=1; j < 32; j++) {
                if ((pos_bits[i] & (1ULL << 2*j)) != 0) { variables.push_back(j+i*32); }
                if ((neg_bits[i] & (1ULL << 2*j)) != 0) { variables.push_back(-j-i*32); }
            }
        }
        return variables;
    }

    std::string to_string() const{
        std::ostringstream oss;
        for (const int& variable : to_variables()) {
            oss << std::to_string(variable);
            oss << ' ';
        }
        return oss.str();
    }

    bool operator==(const bitset &other) const {
        for (int i=0; i<size(); i++) { 
            if (pos_bits[i] != other.pos_bits[i]) return false; 
            if (neg_bits[i] != other.neg_bits[i]) return false; 
        }
        return true;
    }
};

struct visitations {
    // keep track of visited paths
};

struct sat {
    // holds sat problem

    int N;
    std::queue<bit> *unit_queue;
    std::vector<bitset> *clauses;
    std::unordered_map<int, std::vector<std::pair<int, bit>>> *watch_list;

    sat(const int &_N) : N(_N) {
        // allocate onto the heap
        unit_queue = new std::queue<bit>();
        clauses = new std::vector<bitset>();
        watch_list = new std::unordered_map<int, std::vector<std::pair<int, bit>>>();
    }

    ~sat() {
        // delete allocations
        delete unit_queue;
        delete clauses;
        delete watch_list;
    }

    void add_clause(std::unordered_set<int> clause) {
        if (clause.size() == 0) return; // ignore empty
        else if (clause.size() == 1) { // add to unit queue
            bit num(*clause.begin());
            unit_queue->push(num);
            clauses->push_back(bitset(N, clause));
        } else { // add to watchlist 
            bit num1(*clause.begin());
            bit num2(*++clause.begin());
            (*watch_list)[num1.to_variable()].push_back({clauses->size(), num2});
            (*watch_list)[num2.to_variable()].push_back({clauses->size(), num1});
            clauses->push_back(bitset(N, clause));
        }
    }

    int clauses_size() const{ return clauses->size(); }

    int unit_queue_size() const{ return unit_queue->size(); }

    int num_variables() const{ return N; }

    bit pop() {
        bit variable = unit_queue->front();
        unit_queue->pop();
        return variable;
    }

    void emplace(const bit &variable) {
        unit_queue->emplace(variable);
    }

    void clear_queue() {
        while (unit_queue->size() > 0) unit_queue->pop();
    }

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

// debugging -- helper function to verify assignment

bool verify(const bitset &assignment, const std::vector<bitset> &clauses) {
    for (const bitset &clause : clauses) {
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

struct solver {

    // random number generation
    std::mt19937 gen;

    // solver data structures
    sat &sat_data;

    solver(sat &_sat_data) : sat_data(_sat_data) {
        gen = std::mt19937(std::random_device{}());
    }

    // dpll functions

    bool unit_propagation(bitset &assignment) {
        // return conflict 

        // iterate over variables marked for deletion
        while (sat_data.unit_queue_size() > 0) {

            // pop item from queue
            bit bit_num = sat_data.pop();
            if (assignment.test(bit_num.neg())) return false; // early stop

            // assign variable in assignment, cache negative
            assignment.set(bit_num);
            bit bit_num_neg = bit_num.neg();

            // update watched literals
            for (auto pair = (*sat_data.watch_list)[bit_num_neg.to_variable()].begin(); pair != (*sat_data.watch_list)[bit_num_neg.to_variable()].end();) {
                const bitset &bitset_clause = (*sat_data.clauses)[pair->first];
                const bit &bit_other = pair->second;

                // find new watched literal
                bool new_watched_literal_found = false;
                for (int i=0; i<bitset_clause.size(); i++) { // iterate over bitset vector
                    if (new_watched_literal_found) break;

                    // positive bits
                    if (bitset_clause.pos_bits[i] != 0) { // iterate over clause                    
                        for (int j=1; j<32; j++) { // iterate over possible variables in clause
                            if ((bitset_clause.pos_bits[i] & (1ULL<<j*2)) == 0) { // check if variable in clause
                                bit bit_new_num(i*32+j); // construct bit

                                // check bit conditions (not other, negative not assigned)
                                if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num.neg())) {
                                    new_watched_literal_found = true;
                                    (*sat_data.watch_list)[bit_new_num.to_variable()].push_back({pair->first, bit_other}); // add to watchlist

                                    for (std::pair<int, bit> &other_pair : (*sat_data.watch_list)[bit_other.to_variable()]) { // update other
                                        if (other_pair.first == pair->first && other_pair.second == bit_num_neg) {
                                            other_pair.second = bit_new_num;
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }

                    // negative bits
                    if (new_watched_literal_found) break;
                    if (bitset_clause.neg_bits[i] != 0) { // iterate over clause                    
                        for (int j=1; j<32; j++) { // iterate over possible variables in clause
                            if ((bitset_clause.neg_bits[i] & (1ULL<<j*2)) == 0) { // check if variable in clause
                                bit bit_new_num(-i*32-j); // construct bit

                                // check bit conditions (not other, negative not assigned)
                                if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num.neg())) {
                                    new_watched_literal_found = true;
                                    (*sat_data.watch_list)[bit_new_num.to_variable()].push_back({pair->first, bit_other}); // add to watchlist

                                    for (std::pair<int, bit> &other_pair : (*sat_data.watch_list)[bit_other.to_variable()]) { // update other
                                        if (other_pair.first == pair->first && other_pair.second == bit_num_neg) {
                                            other_pair.second = bit_new_num;
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }

                // no new watched literal found
                if (!new_watched_literal_found) {
                    if (assignment.test(bit_other.neg())) return false;
                    else if (!assignment.test_or(bit_other)) sat_data.emplace(bit_other);
                    ++pair;
                } else { // watched literal found
                    pair = (*sat_data.watch_list)[bit_num_neg.to_variable()].erase(pair);
                }
            }
        }

        return true;
    }

    std::string solve(bitset &assignment) {
        // recursive solver

        // visitations... 
        
        // unit propagation
        if (!unit_propagation(assignment)) return "unsat";

        // debugging -- recursive solving
        std::cout << "\tsolving at assignment: " << assignment.to_string() << std::endl;

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

        // debugging -- branching
        std::cout << "\t\tbranching on: " << unassigned_variable << std::endl;

        bitset assignment_true = assignment;
        bit unassigned_bit_true(abs(unassigned_variable));
        sat_data.emplace(unassigned_bit_true);
        std::string result_true = solve(assignment_true);
        if (result_true == "sat") {
            assignment = assignment_true;
            return "sat";
        }

        // debugging -- branching
        std::cout << "\t\tbacktracking and branching on: " << -unassigned_variable << std::endl;
        
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

    static bool test_bits() {
        // tests functionality in bits

        return true;
    }

    static bool test_bitsets() {
        // tests functionality in bitsets

        bitset test_bitset(8);
        test_bitset.set(-1);
        test_bitset.set(-2);
        test_bitset.set(-3);

        if (test_bitset.pos_bits[0] > 0) return false;
        if (test_bitset.neg_bits[0] == 0) return false;

        return true;
    }

    static bool test_sat() {
        // tests functionality in sat
        return true;
    }

    static void run_tests() {
        std::cout << "running unit tests" << std::endl;

        std::cout << "\ttesting bits" << std::endl;
        if (test_bits()) std::cout << "\t\tstatus: passed" << std::endl;
        else std::cout << "\t\tstatus: failed" << std::endl;

        std::cout << "\ttesting bitsets" << std::endl;
        if (test_bitsets()) std::cout << "\t\tstatus: passed" << std::endl;
        else std::cout << "\t\tstatus: failed" << std::endl;

        std::cout << "\ttesting sat" << std::endl;
        if (test_sat()) std::cout << "\t\tstatus: passed" << std::endl;
        else std::cout << "\t\tstatus: failed" << std::endl;
    }

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

    // debugging -- log init variables
    std::cout << "log init variables:" << std::endl;
    std::cout << "\tnumber of variables: " << N << std::endl;
    std::cout << "\tnumber of clauses: " << M << std::endl;

    // // data structures

    // // bitmasks
    //     // position 2i   -> true
    //     // position 2i-1 -> negated
    // // watch list
    //     // variable -> list (clause, other variable)
    
    sat sat_data(N);
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

        // add clause to data structure
        sat_data.add_clause(unorderded_set_clause);
    }

    // debugging -- log reduction
    std::cout << "log reduction:" << std::endl;
    std::cout << "\tmax abs variable: " << max_abs_variable << std::endl;
    std::cout << "\treduction count: " << reduction_count << std::endl;

    // debugging -- stat calculations
    std::cout << "calculated stats" << std::endl;
    std::cout << "\ttotal number of assignments: 2^" << reduction_count << " ~ " << std::pow(2., reduction_count) << std::endl;

    // debugging -- log init information
    std::cout << "log data structure:" << std::endl;
    std::cout << "\tnumber of non unit clauses: " << sat_data.clauses_size() << std::endl;
    std::cout << "\tnumber of unit clauses: " << sat_data.unit_queue_size() << std::endl;
    std::cout << "\twatchlist:\n" << sat_data.watch_list_to_string("\t\t") << std::endl;
    std::cout << "\tclauses:" << std::endl;
    for (const bitset &clause : *sat_data.clauses) {
        std::cout << "\t\tclause: " << clause.to_string() << ", pos bits: ";
        for (int i=0; i< clause.size(); i++) {
            std::cout << clause.pos_bits[i] << ", ";
        }
        std::cout << "neg bits: ";
        for (int i=0; i< clause.size(); i++) {
            std::cout << clause.neg_bits[i] << ", ";
        }
        std::cout << std::endl;
    }

    // debugging -- intermediate time
    auto intermediate = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> intermediate_time = intermediate - start;
    std::cout << "time taken for initialization: " << intermediate_time.count() << std::endl;

    // implementation
    bitset assignment(N);
    solver sat_solver(sat_data);
    std::string result = sat_solver.solve(assignment);
    
    // debugging -- final time
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    const std::chrono::duration<double> since_intermediate_time = end - intermediate;
    std::cout << "total time: " << time.count() << ", since init: " << since_intermediate_time.count() << std::endl;

    // log results
    std::cout << "result: " << result << std::endl;
    if (result == "sat") {
        std::cout << "assignment: " << assignment.to_string() << std::endl;
        std::string verification = verify(assignment, *sat_data.clauses)? "verified" : "failed";
        std::cout << "verification: " << verification << std::endl;
    }

    // default condition
    return 0;

    // notes
        // new changes have broken stuff
            // everything on solveable is positive, which doesn't seem right, but verification works...
            // infeasible returns an assignment, but it's not right...
        // how to debug this stuff...

}
