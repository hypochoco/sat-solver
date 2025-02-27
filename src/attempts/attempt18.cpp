
#include <thread>
#include <iostream>
#include <stack>
#include <queue>
#include <random>

#include "data_structures.h"
#include "watch_list.h"
#include "utilities.h"
#include "vsids.h"

// basics of new multithreading scheme done.. 
// next step is fixing some of the threads
// randomized tasks...

struct solver {

    // variables
    std::mt19937 gen;
    std::stack<std::pair<bit,bitset>> *stack;
    watch_list *wc;
    vsids *decision;
    const std::vector<bitset> &clauses;
    int global_count;
    int base; int iteration_count; int max_iterations;

    solver(const int &_N, const watch_list *_wc, const std::vector<bitset> &_clauses) : clauses(_clauses) {
        gen = std::mt19937(std::random_device{}());
        stack = new std::stack<std::pair<bit,bitset>>();
        wc = new watch_list(*_wc);
        decision = new vsids(_N);
        global_count = 0;
        base = 1e5; iteration_count = 0; max_iterations = 0;
    }

    ~solver() {
        delete stack;
        delete wc;
        delete decision;
    }

    std::optional<bitset> unit_propagation(bitset &assignment) {
        // return conflict 

        // iterate over variables marked for deletion
        while (wc->unit_queue_size() > 0) {

            // pop item from queue
            const bit &bit_num = wc->pop();
            const bit &bit_num_neg = bit_num.neg();

            // assign variable in assignment, cache negative
            assignment.set(bit_num);

            // update watched literals
            std::vector<std::tuple<int,int,bit>> add_list;
            for (auto pair = (*wc->watched_literals)[bit_num_neg.to_variable()].begin(); pair != (*wc->watched_literals)[bit_num_neg.to_variable()].end();) {
                const int bitset_clause_index = pair->first;
                const bitset bitset_clause = clauses[bitset_clause_index];
                const bit bit_other = pair->second;

                // iterate over variables in clause
                bool new_watched_literal_found = false;
                for (const int &num : bitset_clause.to_variables()) { // note: could do randomization here
                    const bit bit_new_num(num);

                    // check bit conditions (not other, negative not assigned)
                    if (!(bit_new_num == bit_other) && !assignment.test(bit_new_num.neg())) {
                        new_watched_literal_found = true;
                        add_list.push_back({bit_new_num.to_variable(),bitset_clause_index,bit_other});

                        for (std::pair<int, bit> &other_pair : (*wc->watched_literals)[bit_other.to_variable()]) { // update other
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
                    if (assignment.test(bit_other.neg())) {
                        // vsids update
                        decision->add(bitset_clause);
                        decision->decay(bitset_clause);

                        // additions
                        for (const std::tuple<int,int,bit> &addition : add_list) {
                            (*wc->watched_literals)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
                        }
                        return bitset_clause;
                    }
                    else if (!assignment.test_or(bit_other)) {
                        if (assignment.test(bit_other.neg())) { return bitset_clause; } // other neg already assigned
                        wc->emplace(bit_other);
                    }
                    ++pair;
                } else { // watched literal found
                    pair = (*wc->watched_literals)[bit_num_neg.to_variable()].erase(pair);
                }
            }
            // additions
            for (const std::tuple<int,int,bit> &addition : add_list) {
                (*wc->watched_literals)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
            }
        }
        return std::nullopt;
    }

    std::string solve(std::atomic<bool> &found, const bit &unassigned_variable, bitset &assignment) {

        // luby restart
        global_count++;
        iteration_count = 0;
        max_iterations = luby(global_count) * base;
        
        // initial iteration
        std::uniform_int_distribution<> dis(0, 1);
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
            const std::optional<bitset> conflict = unit_propagation(assignment);
            if (conflict) { continue; } // skip conflict

            // assign variable
            int unassigned_int = decision->next(assignment); // vsids variable selection
            if (unassigned_int == -1) { found.store(true); return "sat"; } // stopping condition
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
    std::stack<std::pair<bit,bitset>> *stack;
    std::mutex mtx;
    std::vector<std::thread> threads;
    std::atomic<bool> found = false;
    bitset sat_assignment;

    pool(const bitset &default_assignment, std::vector<solver*> _solvers, std::stack<std::pair<bit,bitset>> *_stack) : solvers(_solvers), stack(_stack) {
        std::cout << "started pool with " << solvers.size() << " threads" << std::endl;

        // start tasks
        sat_assignment = bitset(default_assignment);
        for (solver *s : solvers) {
            threads.emplace_back(&pool::task, this, s);
            active_tasks++;
        }
    }

    ~pool() { // wait for all tasks to finish
        for (std::thread &t : threads) if (t.joinable()) t.join(); 
    }

    void task(solver *s) {
        // get assignments from stack

        while (true) {
            bit unassigned_variable;
            bitset assignment = sat_assignment.emptied();
            { // look for assignments, stop if empty
                std::lock_guard<std::mutex> lock(mtx);
                if (stack->size() == 0 || found.load()) { // stop if stack is empty or sat found
                    active_tasks--; // this isn't right... think about this in a sec...
                    return;
                } else {


                    // this shouldn't be a stack... we want to randomly put stuff places...

                    std::pair<bit,bitset> top_pair = stack->top();
                    unassigned_variable = top_pair.first;
                    assignment = top_pair.second;
                    stack->pop();
                }
            }

            // complete assignment
            std::string result = s->solve(found, unassigned_variable, assignment);
            if (result == "sat") {
                std::lock_guard<std::mutex> lock(mtx);
                sat_assignment = assignment;
            } else if (result == "stopped") {
                std::lock_guard<std::mutex> lock(mtx);
                while (s->stack->size() > 0) {
                    stack->emplace(s->stack->top()); s->stack->pop();
                }
            }
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
    int num_solvers = std::thread::hardware_concurrency();

    // preprocess data
    std::tuple<int, watch_list*, std::vector<bitset>*> pp = preprocess(argc, argv);
    int N = std::get<0>(pp);
    watch_list *wc = std::get<1>(pp);
    std::vector<bitset> *clauses = std::get<2>(pp);

    // build solvers
    std::vector<solver*> solvers;
    for (int i=0; i<num_solvers; i++) {
        solvers.emplace_back(new solver(N, wc, *clauses));
    }

    // build problem
    std::stack<std::pair<bit,bitset>> *stack = new std::stack<std::pair<bit,bitset>>();
    stack->push({bit(1),bitset(N)}); stack->push({bit(-1),bitset(N)});

    // instantiate thread pool and execute
    bitset default_assignment(N);
    pool p(default_assignment, solvers, stack);

    // debugging
    while (p.num_active_tasks() > 0) {
        std::cout << "stack size: " << p.stack_size() << ", active tasks: " << p.num_active_tasks() << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    // outputs
    p.await(); 
    std::pair<std::string,bitset> result_pair = p.output();
    std::string result = result_pair.first;
    bitset assignment = result_pair.second;

    std::cout << "result: " << result << ", assignment: " << assignment.to_string() << std::endl;

    if (result == "sat") { // debug verification
        std::string verification = verify(assignment, *clauses)? "passed" : "failed";
        std::cout << "verification: " << verification << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time = end - start;
    std::cout << "total time: " << time.count() << std::endl;

    // clean up
    for (solver *s : solvers) delete s;
    delete stack;
    delete wc;
    delete clauses;

    // default return value
    return 0;
}
