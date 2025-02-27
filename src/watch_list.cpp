
#include <sstream>

#include "decision.h"
#include "watch_list.h"


int watch_list::num_variables() const{ return N; }

bit watch_list::pop() {
    if (unit_queue->size() == 0) throw std::runtime_error("cannot pop with queue size 0");
    bit variable = unit_queue->front();
    unit_queue->pop();
    return variable;
}

void watch_list::sanity_check(const bitset &assignment, const std::string &debug) const{

    // iterate through the watched literals and make sure both are not false

    for (int i=0; i<N; i++) {
        bit num1(i);
        bool num1_neg = assignment.test(num1.neg());
        for (auto it : (*watched_literals)[i]) {
            bit num2 = it.second;

            // both cannot be negative
            if (num1_neg && assignment.test(num2.neg())) throw std::runtime_error("watched literals incorrect: " + debug);
        }
    }
}

bool watch_list::unit_propagation(bitset &assignment, const std::vector<bitset> &clauses, decision &d) {
    // return conflict 

    // iterate over variables marked for deletion
    while (unit_queue_size() > 0) {

        // pop item from queue
        const bit &bit_num = pop();
        const bit &bit_num_neg = bit_num.neg();

        // assign variable in assignment, cache negative
        if (assignment.test(bit_num_neg)) return false;
        assignment.set(bit_num);

        // update watched literals
        std::vector<std::tuple<int,int,bit>> add_list;
        for (auto pair = (*watched_literals)[bit_num_neg.to_variable()].begin(); pair != (*watched_literals)[bit_num_neg.to_variable()].end();) {
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

                    for (std::pair<int, bit> &other_pair : (*watched_literals)[bit_other.to_variable()]) { // update other
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
                    d.add(bitset_clause);
                    d.decay(bitset_clause);

                    // additions
                    for (const std::tuple<int,int,bit> &addition : add_list) {
                        (*watched_literals)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
                    }
                    return false;
                } else if (!assignment.test_or(bit_other)) {
                    emplace(bit_other);
                }
                ++pair;
            } else { // watched literal found
                pair = (*watched_literals)[bit_num_neg.to_variable()].erase(pair);
            }
        }
        // additions
        for (const std::tuple<int,int,bit> &addition : add_list) {
            (*watched_literals)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
        }
    }
    return true;
}

void watch_list::calibrate(const bitset &assignment, const std::vector<bitset> &clauses) {
    // return conflict 

    // iterate over variables marked for deletion
    for (const int &variable : assignment.to_variables()) {
        const bit bit_num(variable);
        const bit bit_num_neg = bit_num.neg();

        // update watched literals
        std::vector<std::tuple<int,int,bit>> add_list;
        for (auto pair = (*watched_literals)[bit_num_neg.to_variable()].begin(); pair != (*watched_literals)[bit_num_neg.to_variable()].end();) {
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

                    for (std::pair<int, bit> &other_pair : (*watched_literals)[bit_other.to_variable()]) { // update other
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
                ++pair;
            } else { // watched literal found
                pair = (*watched_literals)[bit_num_neg.to_variable()].erase(pair);
            }
        }
        // additions
        for (const std::tuple<int,int,bit> &addition : add_list) {
            (*watched_literals)[std::get<0>(addition)].push_back({std::get<1>(addition),std::get<2>(addition)});
        }
    }
}

void watch_list::emplace(const bit &num) { unit_queue->emplace(num); }

void watch_list::clear_queue() {
    while (unit_queue->size() > 0) unit_queue->pop();
}

int watch_list::unit_queue_size() const{ return unit_queue->size(); }

std::string watch_list::watch_list_to_string(const std::string &prefix, const std::vector<bitset> *clauses) const{
    // variable, other, clause
    std::ostringstream oss;
    for (auto &item : *watched_literals) {
        oss << prefix << "variable: " << item.first << std::endl;
        for (auto &pair : item.second) {
            oss << prefix << "\tother: " << pair.second.to_variable() << ", clause: " << (*clauses)[pair.first].to_string() << std::endl;
        }
    }
    return oss.str();
}
