
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <cmath>

// clauses can have both a positve and negative
// split the bitset into two

// helper functions

struct bit {
    // encodes a single variable

    int _i = 0; // location in a bits vector (0->N)
    int _j = 0; // location of bit (2*k, k -> [0,32])
    int _variable = 0; // variable
    uint64_t bits = 0; // bits representation

    bit(const int &variable) { set(variable); }
    bit(const int &i, const int&j, const bool neg) { set(i, j, neg); }

    void set(const int &variable) {
        _variable = variable;
        int value = 2*abs(variable);
        _i = value/64;
        _j = value%64;
        bits |= (1ULL << _j);
        if (variable < 0) bits |= (1ULL << (_j-1));
    }

    void set(const int &i, const int &j, const bool neg) {
        _i = i;
        _j = j;
        _variable = _i*64+_j/2;
        bits |= (1ULL << _j);
        if (neg) {
            bits |= (1ULL << (j-1));
            _variable *= -1;
        }
    }

    bit neg() const{
        bit num_neg = *this;
        num_neg._variable *= -1;
        num_neg.bits ^= 1ULL<<(_j-1);
        return num_neg;
    }

    int to_variable() const {
        return _variable;
    }

    bool operator==(const bit &other) const {
        return _variable == other._variable;
    }
};

struct bitset {
    // encodes multiple variables

    std::vector<uint64_t> bits;

    bitset(const int &N) : bits((2*N + 63) / 64, 0) {}
    bitset(const int &N, const int &variable) : bits((2*N + 63) / 64, 0) {
        set(variable);
    }
    bitset(const int &N, const std::unordered_set<int> &variables) : bits((2*N + 63) / 64, 0) {
        for (const int &variable : variables) set(variable);
    }

    bool test(const int &variable) const{
        int value = 2*abs(variable);
        if ((bits[value/64] & (1ULL << (value%64))) == (1ULL << (value%64))) return true;
        return false;
    }

    bool test_or(const int &variable) const{
        int value = 2*abs(variable);
        if ((bits[value/64] & (1ULL << (value%64))) != 0) return true;
        return false;
    }

    bool test(const bit &variable) const{
        if ((bits[variable._i] & variable.bits) == variable.bits) return true;
        else return false;
    }

    bool test_or(const bit &variable) const{
        if ((bits[variable._i] & variable.bits) != 0) return true;
        else return false;
    }

    void set(const int &variable) {
        int value = 2*abs(variable);
        bits[value/64] |= (1ULL << (value%64));
        if (variable < 0) bits[(value-1)/64] |= (1ULL << ((value-1)%64));
    }

    void set(const bit &variable) {
        bits[variable._i] |= variable.bits;
    }

    int size() const{ return bits.size(); }

    std::vector<int> to_variables() const{
        std::vector<int> values;
        for (int i=0; i<bits.size(); i++) {
            for (int j=1; j < 32; j++) {
                if (bits[i] & (1ULL << 2*abs(j))) {
                    if (bits[i] & (1ULL << (2*abs(j)-1))) values.push_back(-j-i*64);
                    else values.push_back(j+i*64);
                }
            }
        }
        return values;
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
        for (int i=0; i<size(); i++) { if (bits[i] != other.bits[i]) return false; }
        return true;
    }
};

// debugging -- helper functions

bool verify(const bitset &assignment, const std::vector<bitset> &clauses) {
    for (const bitset &clause : clauses) {
        bool satisfied = false;
        for (int i=0; i<assignment.size(); i++) {
            if ((assignment.bits[i] & clause.bits[i]) > 0) {
                satisfied = true;
                break;
            }
        }
        if (!satisfied) return false;
    }
    return true;
}

// dpll functions

bool unit_propagation(
    std::queue<bit> &unit_queue, 
    bitset &assignment, 
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>> &watch_list
) {
    // return conflict 

    // debugging queue
    std::cout << "\tbefore unit prop, queue size: " << unit_queue.size() << std::endl;

    // iterate over variables marked for deletion
    while (unit_queue.size() > 0) {

        // pop item from queue
        bit bit_num = unit_queue.front();
        unit_queue.pop();

        // early stop
        if (assignment.test(bit_num)) return false;

        // // debugging -- watched literal change
        // std::cout << "\t\twithin queue operations!" << std::endl;

        // assign variable, cache negative
        assignment.set(bit_num);
        bit bit_num_neg = bit_num.neg();

        // debugging -- assignment
        std::cout << "\t\tassignment during unit prop: " << assignment.to_string() << std::endl;

        // // debugging -- watched literal change
        // std::cout << "\t\tnegative variable: " << bit_num_neg.to_variable() << std::endl;
        // std::cout << "\t\tsize of watch list: " << watch_list[bit_num.to_variable()].size() << std::endl;
        // std::cout << "\t\tsize of neg watch list: " << watch_list[bit_num_neg.to_variable()].size() << std::endl;
        // std::cout << "\t\twatch list keys: " << std::endl;
        // for (const auto &pair : watch_list) {
        //     std::cout << "\t\t\tkey: " << pair.first << ", equal: " << (pair.first == bit_num_neg.to_variable()) << std::endl;
        // }

        // update watched literals
        for (auto pair = watch_list[bit_num_neg.to_variable()].begin(); pair != watch_list[bit_num_neg.to_variable()].end();) {
            const bitset &bitset_clause = pair->first;
            const bit &bit_other = pair->second;

            // // debugging -- watched literal change
            // std::cout << "\t\trequired to change the watched literal!!" << std::endl;

            // find new watched literal
            bool new_watched_literal_found = false;
            for (int i=0; i<bitset_clause.size(); i++) {
                if (new_watched_literal_found) break;
                if (bitset_clause.bits[i] != 0) { // iterate over clause
                    for (int j=1; j<32; j++) { // iterate over possible variables in clause
                        if ((bitset_clause.bits[i] & (1ULL<<j*2)) == 0) continue; // check if vairable in clause

                        // construct bit
                        bool neg = (bitset_clause.bits[i] & (1ULL<<(j*2-1))) != 0;
                        bit bit_new_num(i, j*2, neg);

                        // // debugging -- bit construction
                        // std::cout << "\t\t\tlooking at: " << bit_new_num.to_variable() << std::endl;
                        
                        // check bit conditions (not other, negative not assigned)
                        if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num.neg())) {
                            new_watched_literal_found = true;
                            watch_list[bit_new_num.to_variable()].push_back({bitset_clause, bit_other});
                            for (std::pair<bitset, bit> &other_pair : watch_list[bit_other.to_variable()]) {
                                if (other_pair.first == bitset_clause && other_pair.second == bit_num_neg) {
                                    other_pair.second = bit_new_num;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                }
            }

            // no new watched literal found
            if (!new_watched_literal_found) {

                if (assignment.test(bit_other.neg())) {

                    // debugging -- failed branch
                    std::cout << "\t\tbranch failed, attempted " << bit_other.to_variable() << " but " << bit_other.neg().to_variable() << " in assignment" << std::endl;
                    std::cout << "\t\t\tcurrent assignment: " << assignment.to_string() << std::endl;
                    
                    return false;
                } else if (!assignment.test_or(bit_other)) {
                    unit_queue.emplace(bit_other);

                    // debugging -- emplaced
                    std::cout << "\t\tclause : " << bitset_clause.to_string() << std::endl;
                    std::cout << "\t\templacing: " << bit_other.to_variable() << std::endl;
                }
                ++pair;
            } else { // watched literal found
                pair = watch_list[bit_num_neg.to_variable()].erase(pair);
            }
        }
    }

    // // debugging -- watchlist
    // std::cout << "\twatch list after unit prop: " << std::endl;
    // for (std::pair<int, std::vector<std::pair<bitset, bit>>> pair : watch_list) {
    //     std::cout << "\t\tvariable: " << pair.first << std::endl;
    //     for (std::pair<bitset, bit> &inner_pair : pair.second) std::cout << "\t\t\tclause: " << inner_pair.first.to_string() << ", other variable: " << inner_pair.second.to_variable() << std::endl;
    // }

    return true;
}

std::string solve(
    const int &num_variables, 
    std::queue<bit> &unit_queue, 
    bitset &assignment, 
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>> &watch_list
) {
    // recursive solver

    // debugging -- init variables
    std::cout << "starting solver" << std::endl;
    std::cout << "\tassignment before unit propagation: " << assignment.to_string() << std::endl;

    // unit propagation
    if (!unit_propagation(unit_queue, assignment, watch_list)) return "unsat";

    // debugging -- after unit propagation
    std::cout << "\tcompleted unit propagation" << std::endl;
    std::cout << "\tassignment after unit propagation: " << assignment.to_string() << std::endl;

    // select unassigned variable
    int unassigned_variable = -1;
    for (int i=1; i<=num_variables; i++) {
        if (!assignment.test_or(i)) {
            unassigned_variable = i;
            break;
        }
    }

    if (unassigned_variable == -1) { // stopping condition
        std::cout << "\tdepleted variables, algorithm complete" << std::endl;
        return "sat";
    }

    // debugging -- selected variable
    std::cout << "\tselected variable: " << unassigned_variable << std::endl;

    // branching

    // note:
        // what about keeping track of branches that we've already seen
        // do we need a copy of the watchlist?

    // debugging -- branching
    std::cout << "\tbranching on " << unassigned_variable << std::endl;

    bitset assignment_true = assignment;
    bit unassigned_bit_true(abs(unassigned_variable));
    unit_queue.push(unassigned_bit_true);
    if (solve(num_variables, unit_queue, assignment_true, watch_list) == "sat") {
        assignment = assignment_true;
        return "sat";
    }

    // debugging -- branching
    std::cout << "\tbranching on " << -unassigned_variable << std::endl;

    // note:
        // backtracking? reassess the watchlist?

    unit_queue = {}; // clear queue form branch
    bitset assignment_false = assignment;
    bit unassigned_bit_false(-abs(unassigned_variable));
    unit_queue.push(unassigned_bit_false);
    if (solve(num_variables, unit_queue, assignment_false, watch_list) == "sat") {
        assignment = assignment_false;
        return "sat";
    }
    return "unsat";
}

// unit tests

int test_unordered_map() {
    // debugging unordered map 
        // hashing doesn't work, switched to to_variable()
        // does not assign the watched literal correctly
    return 0;
}

// main function call

int main(int argc, char* argv[]) {

    // // unit testing
    // test_unordered_map();
    // return 0;

    // plan
        // parse input
        // algorithm
            // inference (step 1, step 2): getting rid of unnecessary variables
                // unit literals
                // pure literals
            // branching (step 3): tree search

        // https://codingnest.com/modern-sat-solvers-fast-neat-and-underused-part-3-of-n/

    // optimized plan
        // bitmaps
        // two watched literals 
        // random restarts

        // fast way of telling if a bit is negative, check n & (n-1) == 0

        // solveable assignment: -1 2 -3 -4 5 6 7 -8 -9 10 -11 12 13 -14 15 16 17 18 -19 -20 21 -22 -23 -24 -25 

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

    // data structures

    // bitmasks
        // position 2i   -> true
        // position 2i-1 -> negated
    // watch list
        // variable -> list (clause, other variable)
    
    bitset assignment(N);
    std::queue<bit> unit_queue;
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>> watch_list;

    // debugging -- verification
    std::vector<bitset> clauses;

    // iterate over input cnf file
    while (std::getline(cnfFile, line)) {
        if (line[0] == 'c') continue; // ignore comments
        if (line[0] == 'p') continue; // ignore stats

        // handle clauses
        std::string value;
        std::stringstream ss(line);
        std::unordered_set<int> unorderded_set_clause; // prevent duplicates, not necessary anymore...         
        while (std::getline(ss, value, ' ')) { // iterate over line
            if (value == "0") continue;
            unorderded_set_clause.insert(std::stoi(value));
        }

        // initialization
            // unit clauses  -> queued for unit propagation
            // unordered set -> bitmask
            // watch list    -> 2 literals per clause
        if (unorderded_set_clause.size() == 0) continue;
        else if (unorderded_set_clause.size() == 1) {
            bit num(*unorderded_set_clause.begin());
            unit_queue.push(num);
        } else {
            bitset bitset_clause(N, unorderded_set_clause);
            bit num1(*unorderded_set_clause.begin());
            bit num2(*++unorderded_set_clause.begin());
            watch_list[*unorderded_set_clause.begin()].push_back({bitset_clause, num2});
            watch_list[*++unorderded_set_clause.begin()].push_back({bitset_clause, num1});
        }

        // debugging -- verification
        clauses.push_back(bitset(N, unorderded_set_clause));
    }

    // implementation here

    // debugging -- initialization
    std::cout << "init stats:" << std::endl;
    std::cout << "\tnum variables: " << stats[0] << ", num clauses: " << stats[1] << std::endl;
    std::cout << "\tunit queue size: " << unit_queue.size() << std::endl;
    std::cout << "\twatch list: " << std::endl;
    for (std::pair<int, std::vector<std::pair<bitset, bit>>> pair : watch_list) {
        std::cout << "\t\tvariable: " << pair.first << std::endl;
        for (std::pair<bitset, bit> &inner_pair : pair.second) std::cout << "\t\t\tclause: " << inner_pair.first.to_string() << ", other variable: " << inner_pair.second.to_variable() << std::endl;
    }

    // solver
    std::string result = solve(N, unit_queue, assignment, watch_list);

    // debugging -- assignment bit
    std::cout << "result: " << result << std::endl;
    std::cout << "final assignment bit: " << assignment.to_string() << std::endl;

    // end timing
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;

    // debugging -- time, verification
    std::cout << "time: " << time.count() << "s" << std::endl;
    if (result == "sat") {
        std::string verified = verify(assignment, clauses)? "verified" : "failed";
        std::cout << "verification: " << verified << std::endl;
    }

    // // output message
    // size_t pos = std::string(argv[1]).find_last_of("/\\");
    // std::string filename = (pos != std::string::npos) ? std::string(argv[1]).substr(pos + 1) : "";
    // std::cout << "file: " << filename << ", time: " << time.count() << "s, result: " << result;
    // if (result == "sat") {
    //     std::cout << ", assignment: " << assignment_to_string(assignment) << std::endl;
    // } else std::cout << std::endl;
    
    // end of implementation
    return 0;
}
