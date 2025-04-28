
#include <iostream>

#include "solver.h"
#include "dpll/dpll.h"

int main(int argc, char* argv[]) {

    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    // solver
    sat_instance instance(argv[1]);
    dpll solver = dpll(instance);
    // drunken_dpll solver = drunken_dpll(instance);
    // solver.solve();
    solver.solve_random();

    // stats
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    std::cout << "time: " << time.count() << ", ";
    if (solver.assignment != "none") {
        std::cout << "sat, ";
        std::string valid = instance.verify(solver.assignment)? "valid" : "invalid";
        std::cout << valid << std::endl;
    } else {
        std::cout << "unsat" << std::endl;
    }

    // default return value
    return 0;
}
