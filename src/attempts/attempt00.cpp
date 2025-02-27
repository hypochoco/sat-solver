
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

void print_clause(const std::vector<int> &clause) {
    for (auto &num : clause) std::cout << num << " ";
    std::cout << std::endl;
}

void print_clauses(const std::vector<std::vector<int>> &clauses) {
    for (auto &clause : clauses) {
        print_clause(clause);
    }
}

std::string collapse_clauses(const std::vector<std::vector<int>> &clauses) {
    std::unordered_set<int> collapsed;
    for (const auto &clause : clauses) {
        for (const auto &num : clause) collapsed.insert(num);
    }
    std::string collapes_clauses = "";
    for (const int &n : collapsed) {collapes_clauses += std::to_string(n) + ' ';}
    return collapes_clauses;
}

std::string order_abs(const std::string &assignment) {
    // order assignment string by abs value

    // early stop
    if (assignment == "") return "";

    // read in assignment
    std::string num;
    std::vector<int> nums;
    std::stringstream ss(assignment);
    while (std::getline(ss, num, ' ')) {
        // std::cout << num << std::endl;
        nums.push_back(std::stoi(num));
    }

    // sort by abs
    std::sort(nums.begin(), nums.end(), [](int a, int b) {return std::abs(a) < std::abs(b);});

    // convert back to string
    std::string ordered_assignment = "";
    for (const int &n : nums) {ordered_assignment += std::to_string(n) + ' ';}
    return ordered_assignment;
}

bool verify(const std::vector<std::vector<int>> &clauses, const std::string &assignment) {
    // return true or false if the given assignment is a truth value

    // std::cout << "attempting to verify.." << std::endl;

    // convert assignment into an arr for access
    std::unordered_set<int> assignment_set;
    std::string var;
    std::stringstream ss(assignment);
    while (std::getline(ss, var, ' ')) {
        assignment_set.insert(std::stoi(var));
    }

    // check assignment per clause
    for (const std::vector<int> &clause : clauses) {
        if (!std::any_of(clause.begin(), clause.end(), [&assignment_set](int num) {
            return assignment_set.find(num) != assignment_set.end();
        })) return false;
    }
    return true;
}

bool remove_unit_literals(std::vector<std::vector<int>> &clauses, std::string &assignment) {
    // remove unit clauses

    // early stop
    const int init_size = clauses.size();
    if (init_size == 0) return false;

    // identify unit clauses for removal
    std::unordered_set<int> remove;
    for (const std::vector<int> &clause : clauses) { 
        if (clause.size() == 1 && remove.count(-clause[0]) == 0) { 
            remove.insert(clause[0]); 

            // add to assignment
            assignment += std::to_string(clause[0]) + ' ';
        } 
    }
    if (remove.size() == 0) return false; // early stop

    // // debugging
    // std::cout << "removing unit literals: ";
    // for (const int &r : remove) std::cout << r << " ";
    // std::cout << std::endl;

    // remove
    clauses.erase(std::remove_if(clauses.begin(), clauses.end(), [&remove](std::vector<int> &clause) {
        // return t/f if x contains an element from remove

        // find clauses for removal
        if (std::any_of(clause.begin(), clause.end(), [&remove](const int &num){
            return remove.find(num) != remove.end();
        })) return true;

        // find variables for removal
        clause.erase(std::remove_if(clause.begin(), clause.end(), [&remove](const int &num) {
            return remove.find(-num) != remove.end();
        }), clause.end());
        return false;
        
    }), clauses.end());

    // // debugging
    // std::cout << "after removal: " << std::endl;
    // print_clauses(clauses);

    return init_size != clauses.size();
}

bool remove_pure_literals(std::vector<std::vector<int>> &clauses, std::string &assignment) {
    // remove pure literals
    const int init_size = clauses.size();

    // early stop
    if (init_size == 0) return false;

    // identify pure literals
    std::unordered_map<int, int> freq_map;
    for (std::vector<int> &clause : clauses) for (int num : clause) freq_map[num]++;
    std::unordered_set<int> remove;
    for (auto &[num, freq] : freq_map) if (freq_map[-num] == 0) {
        remove.insert(num);

        // add to assignment
        assignment += std::to_string(num) + ' ';
    }
    
    // early stop
    if (remove.size() == 0) return false;

    // // debugging
    // std::cout << "removing pure literals: ";
    // for (const int &r : remove) std::cout << r << " ";
    // std::cout << std::endl;

    // remove pure literals clauses
    clauses.erase(std::remove_if(clauses.begin(), clauses.end(), [&remove](std::vector<int> &x) {
        // return t/f if x contains an element from remove

        return std::any_of(x.begin(), x.end(), [&remove](int num) {
            return remove.find(num) != remove.end() || remove.find(-num) != remove.end();
        });
    }), clauses.end());
    return init_size != clauses.size();
}

std::unordered_set<int> possible_branches(std::vector<std::vector<int>> &clauses) {
    // returns all possible nums to branch on

    std::unordered_set<int> nums;
    for (std::vector<int> &clause : clauses) for (int num : clause) nums.insert(num);
    return nums;
}

std::vector<std::vector<int>> branch(std::vector<std::vector<int>> clauses, const int &num) {
    // return modified branch based on assignment

    // branch on num
    clauses.erase(std::remove_if(clauses.begin(), clauses.end(), [&num](std::vector<int> &clause) {
        int i = 0, j = -1;
        for (int &var : clause) {
            if (var == num) {return true;}
            if (var == -num) j = i;
            i++;
        }
        if (j != -1 && j < clause.size()) clause.erase(clause.begin()+j);
        return false;
    }), clauses.end());

    return clauses;
}

bool complete(const std::vector<std::vector<int>> &clauses) {
    // returns if the sat solver is completed or not

    // completions: 
        // 0 or 1 clauses
        // unsat -> any clause is empty

    if (clauses.size() <= 1 || std::any_of(clauses.begin(), clauses.end(), [](const std::vector<int> &clause) {
        return clause.size() == 0;
    })) return true;
    return false;
}

std::string completed_output(const std::vector<std::vector<int>> &clauses) {
    // assuming completed, return variable assignemnts

    std::string assignment = "";
    for (const std::vector<int> &clause : clauses) {
        if (clause.size() == 0) return "none";
        for (const int &num : clause) {
            assignment += std::to_string(num) + ' ';
        }
    }
    return assignment;
}

std::string dfs(std::vector<std::vector<int>> clauses) {
    // dfs

    // check complete
    if (complete(clauses)) return completed_output(clauses);

    // inference
    std::string assignment = "";
    while (remove_unit_literals(clauses, assignment) && remove_pure_literals(clauses, assignment));

    // check complete
    if (complete(clauses)) {
        std::string completed_output_string = completed_output(clauses);
        if (completed_output_string == "none") return completed_output_string;
        else return assignment + completed_output(clauses);
    }

    // dfs
    std::unordered_set<int> possibilities = possible_branches(clauses);
    int num = *possibilities.begin();
    std::string output;
    if ((output = dfs(branch(clauses, num))) != "none") return assignment + std::to_string(num) + ' ' + output;
    if ((output = dfs(branch(clauses, -num))) != "none") return assignment + std::to_string(-num) + ' ' + output;;
    return "none";

    // // debugging
    // auto branch_1 = branch(clauses, num);
    // std::cout << "\tbranching at " << num << std::endl;
    // std::cout << "\t\tcurrent assignment: " << order_abs(assignment) << std::endl;
    // std::cout << "\t\tcollapsed clauses: " << order_abs(collapse_clauses(branch_1)) << std::endl;
    // if ((output = dfs(branch_1)) != "none") return assignment + std::to_string(num) + ' ' + output;

    // auto branch_2 = branch(clauses, -num);
    // std::cout << "\tbacktracking and branching at " << -num << std::endl;
    // std::cout << "\t\tcurrent assignment: " << order_abs(assignment) << std::endl;
    // std::cout << "\t\tcollapsed clauses: " << order_abs(collapse_clauses(branch_2)) << std::endl;
    // if ((output = dfs(branch_2)) != "none") return assignment + std::to_string(-num) + ' ' + output;;
    
    // return "none";
}

int main(int argc, char* argv[]) {

    // plan
        // parse input into a workable format
        // algorithm
            // inference (step 1, step 2): getting rid of unnecessary variables
                // unit literals
                // pure literals
            // branching

    // notes
        // memoization -> comparable structure
        // data structure
            // full, ordered string (elements in clause and whatever else is ordered)

    // optimization notes
        // clauses doesn't have to be a vector, could be an arr of vectors (performance?)
        // freq map could be an arr
        // maybe make some stats on inputs
            // how many elements in a clause? as opposed to the number of clauses, etc.
        // current dfs is not great on memory

        // C200_1806 over time
        // C289_179 over time

    // implementation
    auto start = std::chrono::high_resolution_clock::now();

    // data structures
    int stats[2]; // num variables, num clauses
    std::vector<std::vector<int>> clauses; // 2d vector

    // parse input

    // input file error handling
    if (argc != 2) { std::cerr << "Usage: " << argv[0] << " <input.cnf>" << std::endl; return 1; }
    std::ifstream cnfFile(argv[1]);
    if (!cnfFile) { std::cerr << "Error: Cannot open file " << argv[1] << std::endl; return 1; }

    // iterate over input cnf file
    std::string line;
    while (std::getline(cnfFile, line)) {
        if (line[0] == 'c') continue; // ignore comments
        if (line[0] == 'p') { // initial line
            int count = 0;
            std::string value;
            std::stringstream ss(line);
            while (std::getline(ss, value, ' ')) { // iterate over line
                if (value == "p") continue;
                if (value == "cnf") continue;
                stats[count] = std::stoi(value);
                count++;
            }
        } else {
            std::unordered_set<int> clause_set;
            std::string value;
            std::stringstream ss(line);
            while (std::getline(ss, value, ' ')) { // iterate over line, store into data structure
                if (value == "0") continue;

                clause_set.insert(std::stoi(value));
            }
            clauses.push_back(std::vector<int>(clause_set.begin(), clause_set.end()));
        }
    }

    // implementation here
    std::string result = dfs(clauses);

    // end timing
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;

    // output message
    std::cout << "Execution time: " << time.count() << " seconds, ";
    std::string print_result = (result != "none")? order_abs(result) : "unsat";
    std::cout << "assignment: " << print_result;
    if (result != "none") {
        std::string verification = verify(clauses, result)? "verified sat" : "error--found assignment does not satisfy expression";
        std::cout << ", verification: " << verification << "\n";
    } else {
        std::cout << std::endl;
    }
    
    // end of implementation

    return 0;
}
