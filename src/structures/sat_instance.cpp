
#include <sstream>
#include <fstream>

#include "sat_instance.h"

sat_instance::sat_instance(std::string filepath) {

    // variables
    this->clauses = new std::vector<std::string>();

    // readable file
    std::ifstream cnfFile(filepath);

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
    this->num_variables = stats[0]; // number of variables
    int M = stats[1]; // number of clauses

    // parse file
    while (std::getline(cnfFile, line)) {
        if (line[0] == 'c') continue; // ignore comments
        if (line[0] == 'p') continue; // ignore stats
        std::stringstream ss(line);
        this->clauses->push_back(ss.str());
    }
}

sat_instance::~sat_instance() {
    delete this->clauses;
}

bool sat_instance::verify(std::string &assignment) {

    // read assignment
    std::string value;
    std::stringstream ss(assignment);
    std::vector<int> variables;
    std::vector<bool> values;
    bool flip = true;
    while (ss >> value) {
        if (flip) variables.push_back(std::stoi(value));
        else values.push_back(value == "True");
        flip = !flip;
    }

    // process clauses
    std::vector<std::vector<int>> clauses;
    for (const auto &line : *this->clauses) {
        std::vector<int> clause;
        std::string value;
        std::stringstream ss(line);
        while (ss >> value) {
            if (value == "0") break; 
            int int_value = std::stoi(value);
            clause.push_back(int_value);
        }
        clauses.push_back(clause);
    }

    // check satisfiability
    for (const std::vector<int> clause : clauses) {
        bool satisfied = false;
        for (const int &variable : clause) {
            if (satisfied) break;
            for (int i=0; i<variables.size(); i++) {
                int assigned_variable = variables[i];
                if (!values[i]) assigned_variable *= -1;
                if (variable == assigned_variable){
                    satisfied = true;
                    break;
                }
            }
        }
        if (!satisfied) return false;
    }
    return true;
}
