
#pragma once

#include <unordered_set>
#include <random>
#include <stdexcept>


struct bit {
    int variable = 0;
    bool negative = false;
    int i, j = 0;
    uint64_t bits = 0;

    inline bit() {}
    inline bit(const int &variable) {
        this->variable = variable;
        negative = variable < 0;
        i = abs(variable)/64;
        j = abs(variable)%64;
        bits = 1ULL << j;
    }
    inline int to_variable() const{
        return variable;
    }
    bit neg() const;
    bool operator==(const bit &other) const;
};

struct bitset {
    std::vector<uint64_t> pos_bits;
    std::vector<uint64_t> neg_bits;

    inline bitset() {}
    inline bitset(const int& N) : pos_bits(N/64+1, 0), neg_bits(N/64+1, 0) {}
    inline bitset(const int& N, const std::unordered_set<int> &nums) : pos_bits(N/64+1, 0), neg_bits(N/64+1, 0) {
        for (const int &num : nums) set(num);
    }
    inline bitset(const bitset &other) {
        pos_bits = std::vector<uint64_t>(other.pos_bits);
        neg_bits = std::vector<uint64_t>(other.neg_bits);
    }

    inline void set(const int &num) {
        if (num >= 0) {
            pos_bits[abs(num)/64] |= 1ULL << (abs(num)%64);
        } else {
            neg_bits[abs(num)/64] |= 1ULL << (abs(num)%64);
        }
    }

    void set(const bit &variable);
    void set(const bitset &other);

    void remove(const bit &variable);
    void remove_or(const bit &variable);

    bool test(const bit &bit) const;
    bool test_or(const bit &bit) const;
    bool contains(const bitset &other) const;

    int size() const;
    int num_variables() const;
    bool empty() const;

    bitset neg() const;
    bitset emptied() const;

    std::vector<int> to_variables() const;
    std::string to_string() const;

    bitset operator&(const bitset &other) const;
    bitset operator&=(const bitset &other);
};

template< typename T>
struct random_stack {

    std::vector<T> items;
    std::mt19937 rng;

    inline random_stack() : rng(std::random_device{}()) {}
    inline random_stack(const random_stack &other) : rng(std::random_device{}()) {
        items = other.items;
    }

    inline void emplace(T& x) { items.push_back(x); }
    inline bool empty() const{ return items.empty(); }
    inline T pop_random() {
        if (items.empty()) throw std::runtime_error("random stack is empty");
        std::uniform_int_distribution<size_t> dist(0, items.size() - 1);
        size_t index = dist(rng);
        T result = items[index];
        items[index] = items.back();
        items.pop_back();
        return result;
    }
    inline size_t size() const { return items.size(); }
};