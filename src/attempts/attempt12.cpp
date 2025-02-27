
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

// moving into cdcl and visitations
// resolved vsids, which did a lot...

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
    bitset(const bitset &other) { // copy constructor
        pos_bits = std::vector<uint64_t>(other.pos_bits);
        neg_bits = std::vector<uint64_t>(other.neg_bits);
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

    bitset neg() {
        bitset negated = *this;
        negated.pos_bits = std::vector<uint64_t>(neg_bits);
        negated.neg_bits = std::vector<uint64_t>(pos_bits);
        return negated;
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
    std::unordered_map<int, std::vector<std::pair<int, bit>>> *watch_list;

    sat_problem(const int &_N) : N(_N) {
        // // debugging -- logging sat problem constructor
        // std::cout << "sat problem constructor" << std::endl;

        // heap allocations
        unit_queue = new std::queue<bit>();
        clauses = new std::vector<bitset>();
        watch_list = new std::unordered_map<int, std::vector<std::pair<int, bit>>>();
    }

    ~sat_problem() {
        // // debugging -- logging sat problem destructor
        // std::cout << "sat problem destructor" << std::endl;

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

    for (const bitset &clause : *sat_data.clauses) {
        // check all variables in clause are in the assignment
        std::vector<int> clause_variables = clause.to_variables();
        bool all_variables_assigned = true;
        for (const int &clause_variable : clause_variables) {
            bit bit_clause_variable(clause_variable);
            if (!assignment.test_or(bit_clause_variable)) {
                all_variables_assigned = false;
                break;
            }
        }
        if (!all_variables_assigned) continue;

        // check correct assignment
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
    int c; int j; float c_add; float c_decay;
    std::unordered_map<int, float>* variable_scores; // variable -> value
    std::set<std::pair<float, int>>* sorted_variable_scores; // value -> variable

    vsids(const int &N) {
        // // debugging -- logging vsids constructor
        // std::cout << "vsids constructor" << std::endl;

        // note
            // neither of these two make a large difference
            // having vsids makes a decent difference though

        c = 0;
        // j = 10; c_add = 1; c_decay = 0.95;
        j = 1000; c_add = 1; c_decay = 0.5;
        variable_scores = new std::unordered_map<int, float>();
        sorted_variable_scores = new std::set<std::pair<float, int>>();

        // vairable initialization
        for (int i=0; i<N; i++) {
            (*variable_scores)[i+1] = 0;
            sorted_variable_scores->insert({0, i+1});
        }
    }

    ~vsids() {
        // // debugging -- logging vsids destructor
        // std::cout << "vsids destructor" << std::endl;

        delete variable_scores;
        delete sorted_variable_scores;
    }

    void add(const bitset &clause) {
        // give a conflict clause, update all the variables by c_add
        for (const int &variable : clause.to_variables()) {
            int abs_variable = abs(variable);

            // remove old
            sorted_variable_scores->erase({(*variable_scores)[abs_variable], abs_variable});

            // insert new value
            (*variable_scores)[variable] += c_add;
            sorted_variable_scores->insert({(*variable_scores)[abs_variable], abs_variable});
        }
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
                (*variable_scores)[variable] *= c_decay;
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

    sat_solver(sat_problem &_sat_data, visitations &_visited, const int &N) : sat_data(_sat_data), visited(_visited), decision(N) {

        // // debugging -- logging constructor
        // std::cout << "sat solver constructor" << std::endl;

        // random generation
        gen = std::mt19937(std::random_device{}());

        // luby restarts, default values
        base = 1; iteration_count = 0; max_iterations = 1;
    }

    // ~sat_solver() {
    //     // debugging -- logging destructor
    //     std::cout << "sat solver destructor" << std::endl;
    // }

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
                        return false;
                    }
                    else if (!assignment.test_or(bit_other)) sat_data.emplace(bit_other);
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
        return true;
    }

    std::string solve(bitset &assignment) {
        // recursive solver

        // // debugging -- visitation count
        // std::cout << "visited count: " << visited.count()  << ", visited size: " << visited.visited->size() << std::endl;
        // std::cout << "unit queue size: " << sat_data.unit_queue_size() << std::endl;
        // std::cout << "assignment: " << assignment.to_string() << std::endl;
        
        // luby restarts
        bool luby = true;
        if (luby) {
            if (iteration_count >= max_iterations) return "stopped";
            else iteration_count++;
        }
        
        // visitations
            // this type of pruning does not work
            // we'll have to use backjumping or something else... 
        bool visitations_toggle = true;
        if (visitations_toggle && visited.test(assignment)) {
            return "unsat";
        }
        
        // unit propagation
        if (!unit_propagation(assignment)) return "unsat";

        // // randomized variable selection
        // std::vector<int> unassigned_variables;
        // for (int i=1; i<=sat_data.num_variables(); i++) {
        //     if (!assignment.test_or(i)) { unassigned_variables.push_back(i); }
        // }
        // int unassigned_variable = -1;
        // if (unassigned_variables.size() > 0) {
        //     std::uniform_int_distribution<> dis(0, unassigned_variables.size() - 1);
        //     unassigned_variable = unassigned_variables[dis(gen)];
        // }

        // vsids variable selection
        int unassigned_variable = decision.next(assignment);

        // stopping condition
        if (unassigned_variable == -1) return "sat"; 

        // branching
        bitset assignment_true(assignment);
        bit unassigned_bit_true(abs(unassigned_variable));
        sat_data.emplace(unassigned_bit_true);
        std::string result_true = solve(assignment_true);
        if (result_true == "sat") {
            assignment = assignment_true;
            return "sat";
        } else if (result_true == "stopped") return "stopped";

        bitset visited_assingment_true(assignment);
        visited_assingment_true.set(abs(unassigned_variable));
        if (visitations_toggle) visited.insert(visited_assingment_true);

        sat_data.clear_queue(); // clear queue
        bitset assignment_false(assignment);
        bit unassigned_bit_false(-abs(unassigned_variable));
        sat_data.emplace(unassigned_bit_false);
        std::string result_false = solve(assignment_false);
        if (result_false == "sat") {
            assignment = assignment_false;
            return "sat";
        } else if (result_false == "stopped") return "stopped";

        // bitset visited_assingment_false(assignment);
        // visited_assingment_false.set(-abs(unassigned_variable));
        // if (visitations_toggle) visited.insert(visited_assingment_false);

        // slight improvement
        if (visitations_toggle) visited.insert(assignment);
    
        return "unsat";
    }

    std::string restart(bitset &assignment) {
        
        // notes:
            // what is a visitation? could we not just add in all the clauses?
            // cdcl, adding in learned clauses
            // parallelization
            // vsids

        // C1597_024 -> new target




        // luby restarts
        base = 1e3;
        for (int i=1; i<1e6; i++) {
            std::cout << "luby restart" << std::endl;
            iteration_count = 0;
            max_iterations = luby(i) * base;
            bitset restart_assignment(assignment);
            std::string result = solve(restart_assignment);
            if (result == "sat") { // found result
                assignment = restart_assignment;
                return result;
            } else if (result == "unsat") return result; // found result
        }
        return "stopped";

        // // no luby restarts
        // return solve(assignment);
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

    static void run_tests() {
        std::string test_bit_contains_status = test_bitset_contains()? "passed" : "failed";
        std::cout << "test_bitset_contains: " << test_bit_contains_status << std::endl;
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
    
    // // log results
    // auto end = std::chrono::high_resolution_clock::now();
    // const std::chrono::duration<double> time = end - start;
    // std::cout << "total time: " << time.count() << std::endl;
    // std::cout << "result: " << result << std::endl;
    // if (result == "sat") {
    //     std::cout << "assignment: " << assignment.to_string() << std::endl;

    //     // // partial verification
    //     // std::string partial_verification = partial_verify(assignment, sat_data)? "verified" : "failed";
    //     // std::cout << "partial verification: " << partial_verification << std::endl;

    //     // normal verification
    //     std::string verification = verify(assignment, sat_data)? "verified" : "failed";
    //     std::cout << "verification: " << verification << std::endl;
    // }

    // final output
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    std::string formatted_result = (result == "sat")? "SAT" : "UNSAT";
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
