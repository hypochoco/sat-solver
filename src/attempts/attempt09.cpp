
#include <iostream>
#include <chrono>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <random>

// looks like this is working well, going to remove comments adn test on times...

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
        // std::cout << "adding clause: ";
        // for (const int &num : clause) {
        //     std::cout << num << " ";
        // }
        // std::cout << std::endl;

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

        // now some logic somewhere is not being caught...
        // where is that logic...

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

                    // // debugging -- new watched literal not found
                    // std::cout << "\t\twatched literal for " << bit_num_neg.to_variable() << " not found" << std::endl;

                    if (assignment.test(bit_other.neg())) return false;
                    else if (!assignment.test_or(bit_other)) sat_data.emplace(bit_other);
                    ++pair;
                } else { // watched literal found

                    // // debugging -- new watched literal not found
                    // std::cout << "\t\twatched literal for " << bit_num_neg.to_variable() << " found" << std::endl;
                    // std::cout << "\t\taddlist:" << std::endl;
                    // for (const std::tuple<int,int,bit> &addition : add_list) {
                    //     std::cout << "\t\t\tnew variable: " << std::get<0>(addition) << ", clause index: " <<  std::get<1>(addition) << " , other: " << std::get<2>(addition).to_variable() << std::endl;
                    // }

                    pair = (*sat_data.watch_list)[bit_num_neg.to_variable()].erase(pair);
                }
            }
            // additions
            for (const std::tuple<int,int,bit> &addition : add_list) {
                (*sat_data.watch_list)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
            }
        }

        // // debugging -- watchlist
        // std::cout << "\t\twatchlist after unit prop:" << std::endl;
        // std::cout << sat_data.watch_list_to_string("\t\t\t") << std::endl; 

        return true;
    }

    std::string solve(bitset &assignment) {
        // recursive solver

        // visitations TODO
        
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

    // debugging -- log init variables
    std::cout << "log init variables:" << std::endl;
    std::cout << "\tnumber of variables: " << N << std::endl;
    std::cout << "\tnumber of clauses: " << M << std::endl;

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

    // debugging -- log reduction
    std::cout << "log reduction:" << std::endl;
    std::cout << "\tmax abs variable: " << max_abs_variable << std::endl;
    std::cout << "\treduction count: " << reduction_count << std::endl;

    // debugging -- stat calculations
    std::cout << "calculated stats" << std::endl;
    std::cout << "\ttotal number of assignments: 2^" << reduction_count << " ~ " << std::pow(2., reduction_count) << std::endl;

    // // debugging -- log init information
    // std::cout << "log data structure:" << std::endl;
    // // std::cout << "\tnumber of non unit clauses: " << sat_data.clauses_size() << std::endl;
    // // std::cout << "\tnumber of unit clauses: " << sat_data.unit_queue_size() << std::endl;
    // std::cout << "\twatchlist:\n" << sat_data.watch_list_to_string("\t\t") << std::endl;
    // std::cout << "\tclauses:" << std::endl;
    // for (const bitset &clause : *sat_data.clauses) {
    //     std::cout << "\t\tclause: " << clause.to_string() << ", pos bits: ";
    //     for (int i=0; i< clause.size(); i++) {
    //         std::cout << clause.pos_bits[i] << ", ";
    //     }
    //     std::cout << "neg bits: ";
    //     for (int i=0; i< clause.size(); i++) {
    //         std::cout << clause.neg_bits[i] << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    // debugging -- intermediate time
    auto intermediate = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> intermediate_time = intermediate - start;
    std::cout << "time taken for initialization: " << intermediate_time.count() << std::endl;

    // implementation
    bitset assignment(N);
    sat_solver sat_solver(sat_data);
    std::string result = sat_solver.restart(assignment);
    
    // debugging -- final time
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    const std::chrono::duration<double> since_intermediate_time = end - intermediate;
    std::cout << "total time: " << time.count() << ", since init: " << since_intermediate_time.count() << std::endl;

    // log results
    std::cout << "result: " << result << std::endl;
    if (result == "sat") {
        std::cout << "assignment: " << assignment.to_string() << std::endl;
        std::string verification = verify(assignment, sat_data)? "verified" : "failed";
        std::cout << "verification: " << verification << std::endl;
    }

    // default condition
    return 0;
}
