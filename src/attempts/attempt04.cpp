
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

// this is pretty darn fast, but it looks like its getting stuck in bad search spaces
    // this is apparent as the size of the visitation is pretty small
// adding some randomization will hopefully be beneficial

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
        _variable = _i*32+_j/2;
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

    int c;
    std::vector<uint64_t> pos_bits;
    std::vector<uint64_t> neg_bits;

    bitset(const int &N) : pos_bits((2*N + 63) / 64, 0), neg_bits((2*N + 63) / 64, 0) {}
    bitset(const int &N, const int &variable) : pos_bits((2*N + 63) / 64, 0), neg_bits((2*N + 63) / 64, 0) {
        set(variable);
    }
    bitset(const int &N, const std::unordered_set<int> &variables) : pos_bits((2*N + 63) / 64, 0), neg_bits((2*N + 63) / 64, 0) {
        for (const int &variable : variables) set(variable);
    }

    void set(const int &variable) {
        c++;
        int value = 2*abs(variable);
        if (variable >= 0) {
            pos_bits[value/64] |= (1ULL << (value%64));
        } else {
            neg_bits[value/64] |= (1ULL << (value%64));
            neg_bits[value/64] |= (1ULL << ((value-1)%64));
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
        else return false;
    }

    int size() const{ return pos_bits.size(); }

    int count() const{ return c; }

    std::vector<int> to_variables() const{
        std::vector<int> values;
        for (int i=0; i<pos_bits.size(); i++) {
            for (int j=1; j < 32; j++) {
                if ((pos_bits[i] & (1ULL << 2*j)) != 0) { values.push_back(j+i*32); }
                if ((neg_bits[i] & (1ULL << 2*j)) != 0) { values.push_back(-j-i*32); }
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

    bool contains(const bitset &other) const {
        // this contains other

        for (int i=0; i<size(); i++) {
            if ((pos_bits[i] & other.pos_bits[i]) != other.pos_bits[i]) return false;
            if ((neg_bits[i] & other.neg_bits[i]) != other.neg_bits[i]) return false;
        }
        return true;
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

    int c = 0;
    std::vector<bitset> visited;

    void insert(const bitset &assignment) {
        c++;
        // insert assignment into visited, remove supersets
        for (auto bits = visited.begin(); !(bits == visited.end());) { // remove if needed
            if (bits->contains(assignment)) {
                bits = visited.erase(bits);
            } else ++bits;
        }
        visited.push_back(assignment);
    }

    bool test(const bitset &assignment) {
        for (const bitset &v : visited) { 
            if (assignment.contains(v)) return true;
        }
        return false;
    }

    int size() {
        return visited.size();
    }

    int count() {
        return c;
    }
};

// struct visitations {
//     std::unordered_set<std::string> visited;

//     void insert(const bitset &assignment) {
//         visited.insert(assignment.to_string());
//     }

//     bool test(const bitset &assignment) {
//         if (visited.count(assignment.to_string()) > 0) return true;
//         else return false;
//     }

//     int size() {
//         return visited.size();
//     }
// };

// debugging -- helper functions

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

// dpll functions

bool unit_propagation(
    std::queue<bit> &unit_queue, 
    bitset &assignment, 
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>> &watch_list
) {
    // return conflict 

    // // debugging -- queue
    // std::cout << "\tbefore unit prop, queue size: " << unit_queue.size() << std::endl;

    // iterate over variables marked for deletion
    while (unit_queue.size() > 0) {

        // pop item from queue
        bit bit_num = unit_queue.front();
        unit_queue.pop();

        // early stop
        if (assignment.test(bit_num.neg())) {

            // // debugging -- early stop
            // std::cout << "\t\tbranch failed, attempted " << bit_num.to_variable() << " but " << bit_num.neg().to_variable() << " in assignment: " << assignment.to_string() << std::endl;
            
            return false;
        }

        // // debugging -- watched literal change
        // std::cout << "\t\twithin queue operations!" << std::endl;

        // assign variable, cache negative
        assignment.set(bit_num);
        bit bit_num_neg = bit_num.neg();

        // // debugging -- assignment
        // std::cout << "\t\tassignment during unit prop: " << assignment.to_string() << std::endl;

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
            // // std::cout << "\t\trequired to change the watched literal!!" << std::endl;
            // std::cout << "\t\tupdating watchlist at: " << bit_num_neg.to_variable() << ", other: " << bit_other.to_variable() << ", clause: " << bitset_clause.to_string() << std::endl;

            // find new watched literal
            bool new_watched_literal_found = false;
            for (int i=0; i<bitset_clause.size(); i++) {
                if (new_watched_literal_found) break;

                // positive bits
                if (bitset_clause.pos_bits[i] != 0) { // iterate over clause                    
                    for (int j=1; j<32; j++) { // iterate over possible variables in clause
                        if ((bitset_clause.pos_bits[i] & (1ULL<<j*2)) == 0) continue; // check if vairable in clause

                        // construct bit
                        bit bit_new_num(i, j*2, false);

                        // // debugging -- bit construction
                        // std::cout << "\t\t\tlooking at: " << bit_new_num.to_variable() << std::endl;
                        
                        // check bit conditions (not other, negative not assigned)
                        if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num.neg())) {
                            new_watched_literal_found = true;

                            // // debugging -- new watched literal
                            // std::cout << "\t\t\tnew watched literal found! " << bit_new_num.to_variable() << std::endl;

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

                // negative bits
                if (new_watched_literal_found) break;
                if (bitset_clause.neg_bits[i] != 0) { // iterate over clause                    
                    for (int j=1; j<32; j++) { // iterate over possible variables in clause
                        if ((bitset_clause.neg_bits[i] & (1ULL<<j*2)) == 0) continue; // check if vairable in clause

                        // construct bit
                        bit bit_new_num(i, j*2, true);

                        // // debugging -- bit construction
                        // std::cout << "\t\t\tlooking at: " << bit_new_num.to_variable() << std::endl;
                        
                        // check bit conditions (not other, negative not assigned)
                        if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num.neg())) {
                            new_watched_literal_found = true;

                            // // debugging -- new watched literal
                            // std::cout << "\t\t\tnew watched literal found! " << bit_new_num.to_variable() << std::endl;

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

            // // debugging -- new watched literal
            // if (!new_watched_literal_found) std::cout << "\t\t\tnew watched literal NOT found!" << std::endl;

            // no new watched literal found
            if (!new_watched_literal_found) {

                if (assignment.test(bit_other.neg())) {

                    // // debugging -- failed branch
                    // std::cout << "\t\tbranch failed, attempted " << bit_other.to_variable() << " but " << bit_other.neg().to_variable() << " in assignment" << std::endl;
                    // std::cout << "\t\t\tcurrent assignment: " << assignment.to_string() << std::endl;
                    
                    return false;
                } else if (!assignment.test_or(bit_other)) {
                    unit_queue.emplace(bit_other);

                    // // debugging -- emplaced
                    // std::cout << "\t\t\tnew literal not found in clause : " << bitset_clause.to_string() << ", variable not assigned, emplacing other: " << bit_other.to_variable() << std::endl;
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

std::vector<bit> pure_literal_elimination(
    const int &num_variables,
    const bitset &assignment,
    const std::vector<bitset> &clauses
) {
    // identifies and returns pure literals from unsatisfied clauses based on the given assignment

    // count frequencies
    std::unordered_map<int, int> freq; // variable -> num occurrences
    for (const bitset &clause : clauses) {
        bool satisfied = false;
        for (int i=0; i<assignment.size(); i++) {
            if ((assignment.pos_bits[i] & clause.pos_bits[i]) > 0 || (assignment.neg_bits[i] & clause.neg_bits[i]) > 0) {
                satisfied = true;
                break;
            }
        }
        if (!satisfied) { // clause is unsatisfied
            for (const int variable : clause.to_variables()) {
                freq[variable]++;
            }
        }
    }
    
    // identify pure literals
    std::vector<bit> pure_literals;
    for (int i=1; i<=num_variables; i++) {
        bit num(i);
        if (assignment.test_or(num)) continue; // not in assignment
        if (freq[i] == 0) { pure_literals.push_back(num.neg()); } // pure literal negative
        else if (freq[-i] == 0) { pure_literals.push_back(num); } // pure literal positive
    }
    return pure_literals;
}

std::string solve(
    const int &num_variables, 
    std::queue<bit> &unit_queue, 
    bitset &assignment, 
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>> &watch_list,
    // std::unordered_set<std::string> &visited,
    visitations &visited,
    const std::vector<bitset> &clauses
) {
    // recursive solver
    
    // early stoppage
    if (visited.test(assignment)) return "unsat";

    // // debugging -- init variables
    // std::cout << "starting solver" << std::endl;
    // std::cout << "\tassignment before unit propagation: " << assignment.to_string() << std::endl;

    // unit propagation
    if (!unit_propagation(unit_queue, assignment, watch_list)) return "unsat";

    // unit propagation and pure literals
        // are pure literals worth tracking??
        // pure literals are not really worth tracking...
    // while (true) {
    //     if (!unit_propagation(unit_queue, assignment, watch_list)) return "unsat";
    //     std::vector<bit> pure_literals = pure_literal_elimination(num_variables, assignment, clauses);
    //     if (pure_literals.size() > 0) for (const bit &num : pure_literals) unit_queue.emplace(num);
    //     else break;
    // }

    // // debugging -- after unit propagation
    // std::cout << "\tcompleted unit propagation" << std::endl;
    // std::cout << "\tassignment after unit propagation: " << assignment.to_string() << std::endl;

    // select unassigned variable
    int unassigned_variable = -1;
    for (int i=1; i<=num_variables; i++) {
        if (!assignment.test_or(i)) {
            unassigned_variable = i;
            break;
        }
    }

    if (unassigned_variable == -1) { // stopping condition

        // // debugging -- stoppping condition
        // std::cout << "\tdepleted variables, algorithm complete" << std::endl;

        return "sat";
    }

    // // debugging -- selected variable
    // std::cout << "\tselected variable: " << unassigned_variable << std::endl;

    // branching

    // // debugging -- branching
    // std::cout << "assignment: " << assignment.to_string() << std::endl;

    // note:
        // what about keeping track of branches that we've already seen
        // do we need a copy of the watchlist?

    bitset assignment_true = assignment;
    bit unassigned_bit_true(abs(unassigned_variable));

    // // debugging -- branching
    // std::cout << "\tbranching on " << unassigned_variable << std::endl;
    // std::cout << "true assignment: " << assignment_true.to_string() << std::endl;
    // std::cout << "true assignment: " << assignment_true.to_string() << std::endl;

    unit_queue.push(unassigned_bit_true);
    if (solve(num_variables, unit_queue, assignment_true, watch_list, visited, clauses) == "sat") {
        assignment = assignment_true;
        return "sat";
    }

    visited.insert(assignment_true);
    if (visited.count() % 100 == 0) std::cout << "visited count: " << visited.count() << ", size: " << visited.size() << std::endl;

    // note:
        // backtracking? reassess the watchlist?

    unit_queue = {}; // clear queue form branch
    bitset assignment_false = assignment;
    bit unassigned_bit_false(-abs(unassigned_variable));

    // // debugging -- branching
    // std::cout << "\tbranching on " << -unassigned_variable << std::endl;
    // std::cout << "false assignment: " << assignment_false.to_string() << std::endl;

    unit_queue.push(unassigned_bit_false);
    if (solve(num_variables, unit_queue, assignment_false, watch_list, visited, clauses) == "sat") {
        assignment = assignment_false;
        return "sat";
    }

    visited.insert(assignment_false);
    if (visited.count() % 100 == 0) std::cout << "visited count: " << visited.count() << ", size: " << visited.size() << std::endl;

    return "unsat";
}

// unit tests

int test_bit_and_bitset() {
    // values in the bitset aren't quite right

    bit test_bit_00(536);
    std::cout << "original: " << test_bit_00.to_variable() << std::endl;
    std::cout << "\t_i: " << test_bit_00._i << ", _j: " << test_bit_00._j << std::endl;
    int reconstructed = test_bit_00._i*32+test_bit_00._j/2;
    if ((test_bit_00.bits & test_bit_00.bits-1) != 0) reconstructed *= -1;
    std::cout << "\treconstructed: " << reconstructed << std::endl;
    std::cout << "\tnegative: " << test_bit_00.neg().to_variable() << std::endl;

    bit test_bit_01(-536);
    std::cout << "original: " << test_bit_01.to_variable() << std::endl;
    std::cout << "\t_i: " << test_bit_01._i << ", _j: " << test_bit_01._j << std::endl;
    reconstructed = test_bit_01._i*32+test_bit_01._j/2;
    if ((test_bit_01.bits & test_bit_01.bits-1) != 0) reconstructed *= -1;
    std::cout << "\treconstructed: " << reconstructed << std::endl;
    std::cout << "\tnegative: " << test_bit_01.neg().to_variable() << std::endl;

    return 0;
}

// main function call

int main(int argc, char* argv[]) {

    // // unit testing
    // test_bit_and_bitset();
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
            // luby restarts?

        // visited
            // working well, working fast
            // we are not increasing the size much, meaning that the branches we explore are not that different
            // we need random assignments
            // we need random restarts

    // C175_145
        // conflicts: 220996
        // v -1 -2 -3 -4 -5 -6 7 -8 -9 10 -11 -12 -13 14 -15 -16 -17 18 -19 20 -21 22 23 24 
        // v -25 -26 -27 -28 29 -30 -31 32 33 -34 -35 -36 37 38 -39 -40 -41 42 43 -44 45 
        // v -46 -47 48 -49 50 51 -52 -53 -54 -55 -56 57 58 59 -60 -61 62 -63 -64 65 -66 
        // v 67 68 69 70 -71 72 73 74 75 -76 77 78 79 80 -81 -82 83 -84 85 -86 87 -88 89 
        // v 90 -91 -92 93 94 95 96 -97 -98 -99 -100 -101 102 103 104 -105 -106 -107 108 
        // v -109 -110 -111 112 -113 -114 115 -116 117 118 119 -120 -121 122 -123 -124 
        // v -125 -126 -127 128 -129 -130 -131 132 133 134 -135 -136 -137 138 -139 -140 
        // v -141 -142 -143 -144 145 -146 -147 148 -149 -150 -151 -152 -153 -154 -155 -156 
        // v 157 -158 159 160 -161 162 163 -164 165 166 -167 -168 -169 -170 -171 172 173 
        // v 174 175 0

    // ideas:
        // Variable Selection:
            // VSIDS (Variable State Independent Decaying Sum) to track variable activity
            // Phase saving to remember successful assignments
            // Branching heuristics that focus on variables in recent conflicts
            // Regularly bump scores of variables involved in conflicts

        // Clause Learning and Management:
            // Aggressive clause deletion based on LBD (Literal Block Distance)
            // Keep only high-quality learned clauses
            // Periodically remove less useful learned clauses
            // Use minimization techniques to keep learned clauses small

        // Restarts:
            // Dynamic restart policies based on conflict counts
            // Luby sequence or exponential backoff schemes
            // Local restart policies that preserve progress
            // Track progress metrics to decide when to restart

        // Parallelization:
            // Portfolio-based parallel solvers running different configurations
            // Clause sharing between parallel instances
            // Load balancing across multiple cores
            // Diversification strategies to explore different parts of search space

        // Preprocessing:
            // Variable elimination (resolution)
            // Subsumption checking
            // Binary clause propagation
            // Strong component detection
            // Symmetry detection and breaking

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
    std::vector<bitset> clauses;
    std::queue<bit> unit_queue;
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>> watch_list;
    
    int reduction_count = 1;
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
        while (std::getline(ss, value, ' ')) { // iterate over line
            if (value == "0") continue;
            int int_value = std::stoi(value);

            if (reduction_map[abs(int_value)] == 0) {
                reduction_map[abs(int_value)] = reduction_count;
                reduction_count++;
            }
            if (int_value >= 0) unorderded_set_clause.insert(reduction_map[abs(int_value)]);
            else unorderded_set_clause.insert(-reduction_map[abs(int_value)]);

            // std::cout << "true value: " << int_value << ", abs mapped value: " << reduction_map[abs(int_value)] << std::endl;
        }
        if (unorderded_set_clause.size() == 0) continue; // ignore empty
        else if (unorderded_set_clause.size() == 1) { // if 1, add to unit queue
            bit num(*unorderded_set_clause.begin());
            unit_queue.push(num);
        } else {
            bitset bitset_clause(N, unorderded_set_clause);
            bit num1(*unorderded_set_clause.begin());
            bit num2(*++unorderded_set_clause.begin());
            watch_list[*unorderded_set_clause.begin()].push_back({bitset_clause, num2});
            watch_list[*++unorderded_set_clause.begin()].push_back({bitset_clause, num1});

            // // debugging -- unordered set to bitset
            // std::cout << "num1: " << num1.to_variable() << ", num2: " << num2.to_variable() << std::endl;
            // std::cout << "\tclause: " << bitset_clause.to_string() << std::endl;
            // std::cout << "\t";
            // for (const auto &value : unorderded_set_clause) std::cout << value << ", ";
            // std::cout << std::endl;

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
    // std::unordered_set<std::string> visited;
    visitations visited;
    std::string result = solve(N, unit_queue, assignment, watch_list, visited, clauses);

    // debugging -- assignment bit
    std::cout << "result: " << result << std::endl;
    // std::cout << "final assignment bit: " << assignment.to_string() << std::endl;

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
