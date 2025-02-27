
#include <thread>
#include <iostream>
#include <stack>
#include <queue>
#include <random>
#include <atomic>
#include <mutex>
#include <iomanip>

#include "decision.h"
#include "data_structures.h"
#include "watch_list.h"
#include "utilities.h"


// settings

// // log 17
// int BASE = 1e5;
// int WAIT_TIME = 100;
// int NUM_WORKERS = 12;

// // log 18
// int BASE = 1e3;
// int WAIT_TIME = 100;
// int NUM_WORKERS = 12;

// // log 19, all moms
// int BASE = 1e5;
// int WAIT_TIME = 100;
// int NUM_WORKERS = 12;

// // log 20, half vsids, half moms
// int BASE = 1e5;
// int WAIT_TIME = 100;
// int NUM_WORKERS = 12;

// // log 21, 10 vsids 2 moms
// int BASE = 1e5;
// int WAIT_TIME = 100;
// int NUM_WORKERS = 12;

// log 21, 10 vsids 2 moms
int BASE = 1e3;
int WAIT_TIME = 50;
int NUM_WORKERS = 12;


struct solver {

    // variables
    std::mt19937 gen;
    std::stack<std::pair<bit,bitset>> *stack;
    watch_list *wc;
    decision *d;
    const std::vector<bitset> &clauses;
    int global_count;
    int base; int iteration_count; int max_iterations;

    solver(const int &_N, const watch_list *_wc, const std::vector<bitset> &_clauses, const std::string &decision_type) : clauses(_clauses) {
        gen = std::mt19937(std::random_device{}());
        stack = new std::stack<std::pair<bit,bitset>>();
        wc = new watch_list(*_wc);
        if (decision_type == "vsids") {
            d = new vsids(_N);
        } else if (decision_type == "moms") {
            d = new moms(_N, clauses);
        } else {
            d = new vsids(_N);
        }
        global_count = 0;
        base = BASE; iteration_count = 0; max_iterations = 0;
    }

    ~solver() {
        delete stack; 
        delete wc; 
        delete d;
    }

    std::string solve(std::atomic<bool> &found, const bit &unassigned_variable, bitset &assignment) {

        // luby restart
        global_count++;
        iteration_count = 0;
        max_iterations = luby(global_count) * base;

        // randomize
        d->perturb();
        
        // initial iteration
        std::uniform_int_distribution<> dis(0, 1);
        wc->calibrate(assignment, clauses); // calibration
        stack->push({unassigned_variable,assignment});

        // solve 
        while (stack->size() > 0) {

            // luby restarts
            iteration_count++;
            if (iteration_count > max_iterations) return "stopped";

            // multithreading stop
            if (found.load()) return "unsat";

            // pop assignment from stack, prepare iteration
            std::pair<bit,bitset> pair = stack->top(); stack->pop();
            bit unassigned = pair.first; assignment = pair.second;
            wc->clear_queue(); wc->emplace(unassigned);

            // unit propagation
            const bool conflict = wc->unit_propagation(assignment, clauses, *d);
            if (!conflict) { continue; } // skip conflict

            // assign variable
            int unassigned_int = d->next(assignment); // vsids variable selection
            if (unassigned_int == -1) {
                found.store(true); return "sat"; 
            } // stopping condition
            int polarity = dis(gen) * 2 - 1;
            stack->push({bit(polarity*unassigned_int), assignment});
            stack->push({bit(-polarity*unassigned_int), assignment});
        }
        return "unsat";
    }
};

struct pool {
    
    // variables
    std::atomic<int> active_tasks = 0;
    std::vector<solver*> solvers;
    random_stack<std::pair<bit,bitset>> *stack;
    std::mutex mtx;
    std::vector<std::thread> threads;
    std::atomic<bool> found = false;
    bitset sat_assignment;

    pool(const bitset &default_assignment, std::vector<solver*> _solvers, random_stack<std::pair<bit,bitset>> *_stack) : solvers(_solvers), stack(_stack) {
        sat_assignment = bitset(default_assignment); // start tasks
        for (solver *s : solvers) {
            threads.emplace_back(&pool::task, this, s);
        }
    }

    ~pool() { // wait for all tasks to finish
        for (std::thread &t : threads) if (t.joinable()) t.join(); 
    }

    void task(solver *s) {
        // get assignments from stack

        while (true) {
            if (found.load()) return; // early stop

            bit unassigned_variable;
            bitset assignment = sat_assignment.emptied();

            bool wait = false;
            { // look for assignments
                std::lock_guard<std::mutex> lock(mtx);
                if (stack->size() > 0) { // there are tasks to complete
                    active_tasks++;
                    std::pair<bit,bitset> top_pair = stack->pop_random();
                    unassigned_variable = top_pair.first;
                    assignment = top_pair.second;
                } else if (active_tasks.load() > 0) wait = true; // wait    
                else return; // complete
            }
            if (wait) {
                std::this_thread::sleep_for(std::chrono::milliseconds(WAIT_TIME));
                continue;
            }

            // complete assignment
            std::string result = s->solve(found, unassigned_variable, assignment);
            if (result == "sat") {
                std::lock_guard<std::mutex> lock(mtx);
                sat_assignment = assignment;
            } else if (result == "stopped") {
                std::lock_guard<std::mutex> lock(mtx);
                while (s->stack->size() > 0) {
                    std::pair<bit,bitset> pair = s->stack->top();
                    stack->emplace(pair);
                    s->stack->pop();
                }
            }
            active_tasks--;
        }
    }

    int stack_size() { // thread safe stack size 
        std::lock_guard<std::mutex> lock(mtx);
        return stack->size();
    }

    int num_active_tasks() { // thread safe num active tasks 
        return active_tasks.load();
    }

    void await() {
        for (std::thread &t : threads) if (t.joinable()) t.join();
    }

    std::pair<std::string, bitset> output() { // thread safe output
        std::lock_guard<std::mutex> lock(mtx);
        if (found) {
            return {"sat", sat_assignment};
        } else return {"unsat", sat_assignment};
    }

};

int main(int argc, char* argv[]) {

    // start timer
    auto start = std::chrono::high_resolution_clock::now();

    // settings
    // int num_solvers = std::thread::hardware_concurrency();
    int num_solvers = NUM_WORKERS;

    // preprocess data
    std::tuple<int, watch_list*, std::vector<bitset>*, std::unordered_map<int, int>> pp = preprocess(argc, argv);
    int N = std::get<0>(pp);
    watch_list *wc = std::get<1>(pp);
    std::vector<bitset> *clauses = std::get<2>(pp);
    std::unordered_map<int, int> rr_map = std::get<3>(pp);

    // build solvers
    std::vector<solver*> solvers;
    for (int i=0; i<10; i++) {
        solvers.emplace_back(new solver(N, wc, *clauses, "vsids"));
    }
    for (int i=0; i<2; i++) {
        solvers.emplace_back(new solver(N, wc, *clauses, "moms"));
    }

    // build problem
    vsids decision(N);
    bitset init_assignment(N);
    wc->unit_propagation(init_assignment, *clauses, decision);
    random_stack<std::pair<bit,bitset>> *stack = new random_stack<std::pair<bit,bitset>>();
    int variable = decision.next(init_assignment);
    std::pair<bit,bitset> p1{bit(variable),init_assignment}; std::pair<bit,bitset> p2{bit(-variable),init_assignment};
    stack->emplace(p1); stack->emplace(p2);

    // instantiate thread pool and execute
    bitset default_assignment(N);
    pool p(default_assignment, solvers, stack);

    // // debugging -- log stats
    // while (p.num_active_tasks() > 0) {
    //     std::cout << "stack size: " << p.stack_size() << ", active tasks: " << p.num_active_tasks() << std::endl;
    //     std::this_thread::sleep_for(std::chrono::seconds(1));
    // }

    // outputs
    p.await(); 
    std::pair<std::string,bitset> result_pair = p.output();
    std::string result = result_pair.first;
    bitset assignment = result_pair.second;

    // final output
    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    std::string formatted_result = "";
    if (result == "sat") formatted_result = "SAT";
    else if (result == "unsat") formatted_result = "UNSAT";
    else if (result == "stopped") formatted_result = "STOPPED";
    size_t pos = std::string(argv[1]).find_last_of("/\\");
    std::string filename = (pos != std::string::npos) ? std::string(argv[1]).substr(pos + 1) : "";
    std::cout << "{\"Instance\": \"" << filename << "\", \"Time\": " << std::setprecision(2) << std::fixed << time.count() << ", \"Result\": \"" << formatted_result << "\"";
    if (result == "sat") {
        std::string verification = verify(assignment, *clauses)? "Passed" : "Failed";
        std::cout << ", \"Verification\": \"" << verification << "\"";

        // // reverse reduction map
        // bitset r_assignment = assignment.emptied();
        // for (const int &variable : assignment.to_variables()) {
        //     int r_variable = rr_map[abs(variable)];
        //     r_assignment.set(bit(variable>0? r_variable : -r_variable));
        // }
        // std::stringstream oss;
        // for (const int &variable : r_assignment.to_variables()) {
        //     std::string value = variable > 0? " true " : " false ";
        //     oss << abs(variable) << value;
        // }
        // std::cout << ", \"Solution\": \"" << oss.str() << "\"";
    }
    std::cout << "}" << std::endl;
    
    // clean up
    for (solver *s : solvers) delete s;
    delete stack;
    delete wc;
    delete clauses;

    // default return value
    return 0;
}
