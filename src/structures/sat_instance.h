
#pragma once

#include <iostream>
#include <vector>

struct sat_instance {
    sat_instance(std::string filepath);
    ~sat_instance();
    bool verify(std::string &assignment);

    int num_variables;
    std::vector<std::string>* clauses;
};
