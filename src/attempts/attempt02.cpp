
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

// somewhat working but not really
// branched off of this to rework the bit

// helper functions

struct bit {
    // encodes a single variable

    std::vector<uint64_t> bits;

    bit(const int &N) : bits((2*N + 63) / 64, 0) {};
    bit(const int &N, const int &variable) : bits((2*N + 63) / 64, 0) {
        set(variable);
    }
    bit(const int &N, const int &i, const int&j, const bool neg) : bits((2*N + 63) / 64, 0) {
        set(i, j, neg);
    }

    void set(const int &variable) {
        int value = 2*abs(variable);
        bits[value/64] |= (1ULL << (value%64));
        if (variable < 0) bits[(value-1)/64] |= (1ULL << ((value-1)%64));
    }

    void set(const int &i, const int &j, const bool neg) {
        if (i >= bits.size()) throw std::runtime_error("i>bits.size()");
        bits[i] |= (1ULL << ((j*2)%64));
        if (neg) bits[i] |= (1ULL << (((j*2)-1)%64));
    }

    bit neg() const{
        bit num = *this;
        for (int i=0; i<num.size(); i++) {
            if (num.bits[i] != 0) {
                for (int j=1; j<32; j++) {
                    if (((1 << 2*j) & num.bits[i]) != 0) {
                        num.bits[i] ^= (1ULL<<2*j)-1;
                        return num;
                    }
                }
                return num;
            }
        }
        return num;
    }

    int to_variable() const {
        int value = 0;
        for (int i=0; i<bits.size(); i++) {
            if (bits[i] != 0) {
                double bitslog2 = std::log2(bits[i]); // this is a slow operation
                value = static_cast<int>(bitslog2);
                if (bitslog2 > value) value = -value-i*64;    
                else value += i*64;
                break; // early stop
            }
        }
        return static_cast<int>(value/2);
    }

    int size() const{
        return bits.size();
    }

    bool operator==(const bit &other) const {
        for (int i=0; i<size(); i++) {
            if (bits[i] != 0) {
                return bits[i] == other.bits[i];
            } else return false;
        }
        return false;
    }
};

// struct bithash {
//     std::size_t operator()(const bit &b) const {
//         return std::hash<int>()(b.to_variable());
//     }
// };

struct bithash {
    std::size_t operator()(const bit &b) const {
        std::size_t hash_value = 0;
        for (auto bit_val : b.bits) {
            hash_value ^= std::hash<uint64_t>()(bit_val) + 0x9e3779b9 + (hash_value << 6) + (hash_value >> 2);
        }
        return hash_value;
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
        // NOTE: does not check positive/negative

        int value = 2*abs(variable);
        if (bits[value/64] & (1ULL << (value%64))) return true;
        return false;
    }

    bool test(const bit &variable) const{
        // NOTE: does not check positive/negative

        for (int i=0; i<variable.size(); i++) {
            if (variable.bits[i] != 0) {
                if ((variable.bits[i] & bits[i]) != 0) return true;
                else return false;
            }
        }
        return false;
    }

    bool test_neg(const bit &variable) const{
        // NOTE: only checks negative

        for (int i=0; i<variable.size(); i++) {
            if (variable.bits[i] != 0) {
                uint64_t value = variable.bits[i] & bits[i];
                if (variable.bits[i] == value) return true;
                else return false;
            }
        }
        return false;
    }

    void set(const int &variable) {
        int value = 2*abs(variable);
        bits[value/64] |= (1ULL << (value%64));
        if (variable < 0) bits[(value-1)/64] |= (1ULL << ((value-1)%64));
    }

    void set(const bit &variable) {
        for (int i=0; i<variable.size(); i++) {
            if (variable.bits[i] != 0) {
                bits[i] |= variable.bits[i];
                return;
            }
        }
    }

    int size() const{
        return bits.size();
    }

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
        for (int i=0; i<bits.size(); i++) {
            for (int j=1; j < 32; j++) {
                if (bits[i] & (1ULL << 2*abs(j))) {
                    if (bits[i] & (1ULL << (2*abs(j)-1))) {
                        oss << std::to_string(-j-i*64);
                        oss << ' ';
                    } else {
                        oss << std::to_string(j+i*64);
                        oss << ' ';
                    }
                }
            }
        }
        return oss.str();
    }

    bool operator==(const bitset &other) const {
        for (int i=0; i<size(); i++) {
            if (bits[i] != other.bits[i]) return false;
        }
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
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>, bithash> &watch_list
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
        if (assignment.test_neg(bit_num)) return false;

        // debugging -- watched literal change
        std::cout << "\t\twithin queue operations!" << std::endl;

        // assign variable
        assignment.set(bit_num);

        // update watched literals
        bit bit_num_neg = bit_num.neg();

        // debugging -- watched literal change
        std::cout << "\t\tnegative variable: " << bit_num_neg.to_variable() << std::endl;
        std::cout << "\t\tsize of watch list: " << watch_list[bit_num.to_variable()].size() << std::endl;
        std::cout << "\t\tsize of neg watch list: " << watch_list[bit_num_neg.to_variable()].size() << std::endl;
        std::cout << "\t\twatch list keys: " << std::endl;
        for (const auto &pair : watch_list) {
            std::cout << "\t\t\tkey: " << pair.first << ", equal: " << (pair.first == bit_num_neg.to_variable()) << std::endl;
        }

        for (auto pair = watch_list[bit_num_neg.to_variable()].begin(); pair != watch_list[bit_num_neg.to_variable()].end();) {
            const bitset &bitset_clause = pair->first;
            const bit &bit_other = pair->second;

            // debugging -- watched literal change
            std::cout << "\t\trequired to change the watched literal!!" << std::endl;

            // find new watched literal
            bool new_watched_literal_found = false;
            for (int i=0; i<bitset_clause.size(); i++) {
                if (new_watched_literal_found) break;
                if (bitset_clause.bits[i] != 0) { // iterate over clause
                    for (int j=1; j<32; j++) {

                        if ((bitset_clause.bits[i] & (1ULL<<j*2)) == 0) continue;
                        bool neg = (bitset_clause.bits[i] & (1ULL<<(j*2-1))) != 0;
                        bit bit_new_num(bit_num.size(), i, j, neg);

                        std::cout << "\t\t\tlooking at: " << bit_new_num.to_variable() << std::endl;
                        std::cout << "\t\t\tderping at: " << i*64+j << std::endl;
                        std::cout << "\t\t\tderping at: " << neg << std::endl;

                        if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num)) {
                            new_watched_literal_found = true;
                            watch_list[bit_new_num.to_variable()].push_back({bitset_clause, bit_other});
                            for (std::pair<bitset, bit> &other_pair : watch_list[bit_other.to_variable()]) {
                                if (other_pair.first == bitset_clause && other_pair.second == bit_num_neg) {
                                    other_pair.second = bit_new_num;
                                    break;
                                }
                            }
                            break;
                        } else {
                            bit bit_new_num_neg = bit_new_num.neg();
                            if (!(bit_new_num_neg == bit_other) && !assignment.test(bit_new_num_neg)) {
                                new_watched_literal_found = true;
                                watch_list[bit_new_num_neg.to_variable()].push_back({bitset_clause, bit_other});
                                for (std::pair<bitset, bit> &other_pair : watch_list[bit_other.to_variable()]) {
                                    if (other_pair.first == bitset_clause && other_pair.second == bit_num_neg) {
                                        other_pair.second = bit_new_num_neg;
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
                if (assignment.test_neg(bit_other.neg())) {
                    std::cout << "branch failed: " << bit_other.neg().to_variable() << "in assignment" << std::endl;
                    return false;
                } else if (!assignment.test(bit_other)) unit_queue.emplace(bit_other); 
            }

            // next new_watched_literal_found
            if (new_watched_literal_found) {
                pair = watch_list[bit_num_neg.to_variable()].erase(pair);
            } else {
                ++pair;
            }
        }
    }

    // debugging -- watchlist
    std::cout << "\twatch list after unit prop: " << std::endl;
    for (std::pair<int, std::vector<std::pair<bitset, bit>>> pair : watch_list) {
        std::cout << "\t\tvariable: " << pair.first << std::endl;
        for (std::pair<bitset, bit> &inner_pair : pair.second) std::cout << "\t\t\tclause: " << inner_pair.first.to_string() << ", other variable: " << inner_pair.second.to_variable() << std::endl;
    }

    return true;
}

std::string solve(
    const int &num_variables, 
    std::queue<bit> &unit_queue, 
    bitset &assignment, 
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>, bithash> &watch_list
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
        if (!assignment.test(i)) {
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

    // debugging -- branching
    std::cout << "\tbranching on " << unassigned_variable << std::endl;

    bitset assignment_true = assignment;
    bit unassigned_bit_true(num_variables, unassigned_variable);
    unit_queue.push(unassigned_bit_true); // do we need a copy of this?
    if (solve(num_variables, unit_queue, assignment_true, watch_list) == "sat") {
        assignment = assignment_true;
        return "sat";
    }

    // debugging -- branching
    std::cout << "\tbranching on " << -unassigned_variable << std::endl;

    bitset assignment_false = assignment;
    bit unassigned_bit_false(num_variables, -unassigned_variable);
    unit_queue.push(unassigned_bit_false);
    if (solve(num_variables, unit_queue, assignment_false, watch_list) == "sat") {
        assignment = assignment_false;
        return "sat";
    }
    return "unsat";
}

int test_unordered_map() {
    // debugging unordered map 
        // hashing doesn't work, switched to to_variable()
        // does not assign the watched literal correctly
    return 0;
}

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

        // unsat is working? but sat says unsat
        // go on a large debugging spree?
            // revisit the logic bitset stuff

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
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>, bithash> watch_list;

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
            bit num(N, *unorderded_set_clause.begin());
            unit_queue.push(num);
        } else {
            bitset bitset_clause(N, unorderded_set_clause);
            bit num1(N, *unorderded_set_clause.begin());
            bit num2(N, *++unorderded_set_clause.begin());
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
