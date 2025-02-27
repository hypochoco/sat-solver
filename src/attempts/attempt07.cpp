
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
#include <random>

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
        return false;
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

// struct visitations {

//     int c = 0;
//     std::vector<bitset> visited;

//     void insert(const bitset &assignment) {
//         c++;
//         // insert assignment into visited, remove supersets
//         for (auto bits = visited.begin(); !(bits == visited.end());) { // remove if needed
//             if (bits->contains(assignment)) {
//                 bits = visited.erase(bits);
//             } else ++bits;
//         }
//         visited.push_back(assignment);
//     }

//     bool test(const bitset &assignment) {
//         for (const bitset &v : visited) { 
//             if (assignment.contains(v)) return true;
//         }
//         return false;
//     }

//     int size() {
//         return visited.size();
//     }

//     int count() {
//         return c;
//     }
// };

struct visitations { // testing visitations
    std::unordered_set<std::string> visited;
    void insert(const bitset &assignment) {}
    bool test(const bitset &assignment) { return false; }
    int size() { return visited.size(); }
    int count() { return size(); }
};

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

// luby restarts

int luby(int k) {
    // starts at 1

    int power = 1;
    while (power <= k + 1) power *= 2;
    power /= 2;
    if (power - 1 == k) return power / 2;
    return luby(k - (power - 1));
}

// sat solver class

struct sat_solver {

    // notes
        // pure literals elimination could be made faster with a freq map
        // randomization could be added to unit propogation

    // random number generation
    std::mt19937 gen;

    // solver data structures
    int &num_variables;
    std::queue<bit> &unit_queue;
    visitations &visited;
    std::unordered_map<int, std::vector<std::pair<bitset, bit>>> &watch_list;
    std::vector<bitset> &clauses;

    // luby restarts
    int base;
    int iteration_count;
    int max_iterations;

    sat_solver(
        int &_num_variables,
        std::queue<bit> &_unit_queue,
        visitations &_visited,
        std::unordered_map<int, std::vector<std::pair<bitset, bit>>> &_watch_list,
        std::vector<bitset> &_clauses
    ) : num_variables(_num_variables),
        unit_queue(_unit_queue),
        visited(_visited),
        watch_list(_watch_list),
        clauses(_clauses) 
    {
        // random generation
        gen = std::mt19937(std::random_device{}());

        // luby restarts, default values
        base = 1; iteration_count = 0; max_iterations = 1;
    }

    // dpll functions

    bool unit_propagation(bitset &assignment) {
        // return conflict 

        // iterate over variables marked for deletion
        while (unit_queue.size() > 0) {

            // pop item from queue
            bit bit_num = unit_queue.front();
            unit_queue.pop();

            // early stop
            if (assignment.test(bit_num.neg())) return false;

            // assign variable, cache negative
            assignment.set(bit_num);
            bit bit_num_neg = bit_num.neg();

            // update watched literals
            for (auto pair = watch_list[bit_num_neg.to_variable()].begin(); pair != watch_list[bit_num_neg.to_variable()].end();) {
                const bitset &bitset_clause = pair->first;
                const bit &bit_other = pair->second;

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

                    // negative bits
                    if (new_watched_literal_found) break;
                    if (bitset_clause.neg_bits[i] != 0) { // iterate over clause                    
                        for (int j=1; j<32; j++) { // iterate over possible variables in clause
                            if ((bitset_clause.neg_bits[i] & (1ULL<<j*2)) == 0) continue; // check if vairable in clause

                            // construct bit
                            bit bit_new_num(i, j*2, true);

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
                    if (assignment.test(bit_other.neg())) return false;
                    else if (!assignment.test_or(bit_other)) unit_queue.emplace(bit_other);
                    ++pair;
                } else { // watched literal found
                    pair = watch_list[bit_num_neg.to_variable()].erase(pair);
                }
            }
        }

        return true;
    }

    std::vector<bit> pure_literal_elimination(const bitset &assignment) {
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

    std::string solve(bitset &assignment, bool luby) {
        // recursive solver

        // luby restart early stop
        if (luby) {
            if (iteration_count >= max_iterations) return "stopped";
            else iteration_count++;
        }
        
        // early stoppage
        if (visited.test(assignment)) return "unsat";

        // unit propagation
        if (!unit_propagation(assignment)) return "unsat";

        // unit propagation and pure literals, pure literals are currently not worth tracking
        // while (true) {
        //     if (!unit_propagation(unit_queue, assignment, watch_list)) return "unsat";
        //     std::vector<bit> pure_literals = pure_literal_elimination(num_variables, assignment, clauses);
        //     if (pure_literals.size() > 0) for (const bit &num : pure_literals) unit_queue.emplace(num);
        //     else break;
        // }

        // randomized variable selection
        std::vector<int> unassigned_variables;
        for (int i=1; i<=num_variables; i++) {
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
        std::cout << "assignment: " << assignment.to_string() << std::endl;



        bitset assignment_true = assignment;
        bit unassigned_bit_true(abs(unassigned_variable));
        unit_queue.push(unassigned_bit_true);
        std::string result_true = solve(assignment_true, luby);
        if (result_true == "sat") {
            assignment = assignment_true;
            return "sat";
        } else if (result_true == "stopped") return "stopped";
        visited.insert(assignment_true);

        // // debugging -- visitations
        // if (visited.count() % 10000 == 0) std::cout << "visited count: " << visited.count() << ", size: " << visited.size() << std::endl;

        unit_queue = {}; // clear queue form branch
        bitset assignment_false = assignment;
        bit unassigned_bit_false(-abs(unassigned_variable));
        unit_queue.push(unassigned_bit_false);
        std::string result_false = solve(assignment_false, luby);
        if (result_false == "sat") {
            assignment = assignment_false;
            return "sat";
        } else if (result_false == "stopped") return "stopped";
        visited.insert(assignment_false);

        // // debugging -- visitations
        // if (visited.count() % 10000 == 0) std::cout << "visited count: " << visited.count() << ", size: " << visited.size() << std::endl;

        return "unsat";
    }

    std::string restart(bitset &assignment) {

        // stats: U50_1065_045
            // base: 1e2, visited count: 4165100, size: 16120, time: exceeded
            // base: 1e3, visited count: 2231400, size: 1947,  time: 41.512s
            // base: 1e4, visited count: 3617600, size: 488,   time: 52.5583s
            // none:      visited count: 276200,  size: 10,    time: 6.83314s

        // stats: U75_1597_024
            // none: visited count: 2885300, size: 14, time: exceeded

        // // luby restarts
        // base = 1e3;
        // for (int i=1; i<1e6; i++) {
        //     std::cout << "luby restart" << std::endl;
        //     iteration_count = 0;
        //     max_iterations = luby(i) * base;
        //     bitset restart_assignment(assignment.size());
        //     std::string result = solve(restart_assignment, true);
        //     if (result == "sat") {
        //         assignment = restart_assignment;
        //         return result;
        //     } else if (result == "unsat") return result;
        // }
        // return "stopped";

        // no luby restarts
        return solve(assignment, false);
    }

};

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

int test_luby_sequence() {
    // test correctness of luby sequence

    std::cout << "starting luby sequence test" << std::endl;

    std::cout << "luby sequence: ";
    for (int i=0; i<10; i++) {
        std::cout << luby(i) << " ";
    }
    std::cout << std::endl;

    return 0;
}

// main function call

int main(int argc, char* argv[]) {

    // // unit testing
    // test_bit_and_bitset();
    // test_luby_sequence();
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
        // random restarts
            // luby restarts?
            // luby restarts cause a stack overflow... better memory management required...
        // after running out of time, declare unsat?
        // multithreading?

    // C181_3151 -> returned as UNSAT, should be SAT
        // v -1 -2 -3 -4 -5 6 -7 -8 -9 -10 -11 -12 -13 -14 15 -16 -17 -18 -19 -20 -21 -22 
        // v -23 -24 -25 -26 27 -28 -29 -30 -31 -32 -33 34 -35 -36 -37 -38 -39 -40 -41 -42 
        // v -43 -44 -45 -46 -47 48 -49 -50 -51 -52 53 -54 -55 -56 -57 -58 -59 -60 -61 -62 
        // v -63 -64 -65 -66 -67 -68 69 -70 -71 72 -73 -74 -75 -76 -77 -78 -79 -80 -81 -82 
        // v -83 -84 -85 -86 -87 -88 -89 90 91 -92 -93 -94 -95 -96 -97 -98 -99 -100 101 
        // v -102 -103 -104 -105 -106 -107 -108 -109 -110 111 -112 -113 -114 -115 -116 
        // v -117 -118 -119 -120 121 -122 -123 -124 -125 -126 -127 -128 -129 -130 131 -132 
        // v -133 -134 -135 -136 -137 -138 -139 -140 141 -142 -143 -144 -145 -146 -147 
        // v -148 -149 -150 151 -152 -153 -154 -155 -156 -157 -158 -159 -160 161 -162 -163 
        // v -164 -165 -166 -167 -168 -169 -170 171 -172 -173 -174 -175 -176 -177 -178 
        // v -179 -180 181 0

    // C459_4675 -> returned as UNSAT, should be SAT

    // how to go about testing this...
    // maybe some better visualizations or information on what was and wasn't tested
    // what about conditions where the number of variables doesn't match... 
    // having some statistics about that would be helpful...




    // truth values:

    // C175_145 -> SAT
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
    std::queue<bit> unit_queue;
    // std::vector<bitset> clauses;
    // std::unordered_map<int, std::vector<std::pair<bitset, bit>>> watch_list;
    std::unique_ptr<std::vector<bitset>> clauses = std::make_unique<std::vector<bitset>>();
    std::unique_ptr<std::unordered_map<int, std::vector<std::pair<bitset, bit>>>> watch_list = std::make_unique<std::unordered_map<int, std::vector<std::pair<bitset, bit>>>>();

    // note
        // make a struct for this representation?
        // the watchlist could point to the clauses list... 

    // std::unordered_set<std::string> visited;
    visitations visited;
    
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
        }
        if (unorderded_set_clause.size() == 0) continue; // ignore empty
        else if (unorderded_set_clause.size() == 1) { // if 1, add to unit queue
            bit num(*unorderded_set_clause.begin());
            unit_queue.push(num);
        } else {
            bitset bitset_clause(N, unorderded_set_clause);
            bit num1(*unorderded_set_clause.begin());
            bit num2(*++unorderded_set_clause.begin());
            (*watch_list)[*unorderded_set_clause.begin()].push_back({bitset_clause, num2});
            (*watch_list)[*++unorderded_set_clause.begin()].push_back({bitset_clause, num1});
        }
        clauses->push_back(bitset(N, unorderded_set_clause));
    }

    // implementation here
    sat_solver sat(N, unit_queue, visited, *watch_list, *clauses);
    // std::string result = sat.solve(assignment);
    std::string result = sat.restart(assignment);

    // end timing
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;

    // // debugging -- result, time, verification
    // std::cout << "result: " << result << std::endl;
    // std::cout << "assignment: " << assignment.to_string() << std::endl;
    // std::cout << "time: " << time.count() << "s" << std::endl;
    // if (result == "sat") {
    //     std::string verified = verify(assignment, *clauses)? "verified" : "failed";
    //     std::cout << "verification: " << verified << std::endl;
    // }

    // output message
    std::string formatted_result = (result == "sat")? "SAT" : "UNSAT";
    size_t pos = std::string(argv[1]).find_last_of("/\\");
    std::string filename = (pos != std::string::npos) ? std::string(argv[1]).substr(pos + 1) : "";
    std::cout << "{\"Instance\": \"" << filename << "\", \"Time\": " << std::setprecision(2) << std::fixed << time.count() << ", \"Result\": \"" << formatted_result << "\"}" << std::endl;
    
    // end of implementation
    return 0;
}
