
#include <vector>
#include <unordered_set>
#include <sstream>

#include "data_structures.h"

// bit

bit bit::neg() const{
    bit neg_bit = *this;
    neg_bit.variable *= -1;
    neg_bit.negative = !negative;
    return neg_bit;
}

bool bit::operator==(const bit &other) const {
    return 
    variable == other.variable;
}


// bitset

bitset bitset::emptied() const{
    bitset empty_copy = *this;
    for (uint64_t &bit : empty_copy.pos_bits) { bit = 0; }
    for (uint64_t &bit : empty_copy.neg_bits) { bit = 0; }
    return empty_copy;
}

void bitset::set(const bit &variable) {
    if (!variable.negative) { // positive
        pos_bits[variable.i] |= variable.bits;
    } else {
        neg_bits[variable.i] |= variable.bits;
    }
}

void bitset::set(const bitset &other) {
    for (int i=0; i<std::min(size(), other.size()); i++) {
        pos_bits[i] |= other.pos_bits[i];
        neg_bits[i] |= other.neg_bits[i];
    }
}

void bitset::remove(const bit &variable) {
    if (!variable.negative) { // positive
        if ((pos_bits[variable.i] & variable.bits) > 0) pos_bits[variable.i] ^= variable.bits;
    } else {
        if ((neg_bits[variable.i] & variable.bits) > 0) neg_bits[variable.i] ^= variable.bits;
    }
}

void bitset::remove_or(const bit &variable) {
    if ((pos_bits[variable.i] & variable.bits) > 0) pos_bits[variable.i] ^= variable.bits;
    if ((neg_bits[variable.i] & variable.bits) > 0) neg_bits[variable.i] ^= variable.bits;
}

bool bitset::test(const bit &bit) const{
    if (!bit.negative) { // positive
        return (pos_bits[bit.i] & bit.bits) > 0;
    } else { // negative
        return (neg_bits[bit.i] & bit.bits) > 0;
    }
}

bool bitset::test_or(const bit &bit) const{
    return (pos_bits[bit.i] & bit.bits) > 0 || (neg_bits[bit.i] & bit.bits) > 0;
}

bool bitset::contains(const bitset &other) const{
    // returns if this bitset contains the other
    if (other.pos_bits.size() > pos_bits.size()) return false; // early stop
    if (other.neg_bits.size() > neg_bits.size()) return false; // early stop

    for (int i=0; i<pos_bits.size(); i++){ // pos bits
        // other.pos_bits cannot have a bit that pos_bits does not
        if ((pos_bits[i] & other.pos_bits[i]) != other.pos_bits[i]) return false;
    }
    for (int i=0; i<neg_bits.size(); i++){ // neg bits
        // other.pos_bits cannot have a bit that pos_bits does not
        if ((neg_bits[i] & other.neg_bits[i]) != other.neg_bits[i]) return false;
    }

    return true;
}

int bitset::size() const{ return std::max(pos_bits.size(), neg_bits.size()); }

int bitset::num_variables() const{
    int count = 0;
    for (int i=0; i<std::max(pos_bits.size(), neg_bits.size()); i++) {
        if (pos_bits[i] == 0 && neg_bits[i] == 0) continue;
        for (int j=0; j<64; j++) {
            if (pos_bits[i] & (1ULL << j)) count++;
            if (neg_bits[i] & (1ULL << j)) count++;
        }
    }
    return count;
}

bitset bitset::neg() const{
    bitset negated = *this;
    negated.pos_bits = std::vector<uint64_t>(neg_bits);
    negated.neg_bits = std::vector<uint64_t>(pos_bits);
    return negated;
}

bool bitset::empty() const{
    for (int i=0; i<std::max(pos_bits.size(), neg_bits.size()); i++) {
        if (pos_bits[i] > 0 || neg_bits[i] > 0) return false;
    }
    return true;
}

std::vector<int> bitset::to_variables() const{
    std::vector<int> variables;
    for (int i=0; i<std::max(pos_bits.size(), neg_bits.size()); i++) {
        if (pos_bits[i] == 0 && neg_bits[i] == 0) continue;
        for (int j=0; j<64; j++) {
            if (pos_bits[i] & (1ULL << j)) variables.push_back(i*64+j);
            if (neg_bits[i] & (1ULL << j)) variables.push_back(-i*64-j);
        }
    }
    return variables;
}

std::string bitset::to_string() const{
    std::stringstream oss;
    for (int i=0; i<std::max(pos_bits.size(), neg_bits.size()); i++) {
        if (pos_bits[i] == 0 && neg_bits[i] == 0) continue;
        for (int j=0; j<64; j++) {
            if (pos_bits[i] & (1ULL << j)) oss << i*64+j << ' ';
            if (neg_bits[i] & (1ULL << j)) oss << -i*64-j << ' ';
        }
    }
    return oss.str();
}

bitset bitset::operator&(const bitset &other) const{
    bitset copy = *this;
    for (int i=0; i<other.size(); i++) {
        copy.pos_bits[i] &= other.pos_bits[i];
        copy.neg_bits[i] &= other.neg_bits[i];
    }
    return copy;
}

bitset bitset::operator&=(const bitset &other) {
    for (int i=0; i<other.size(); i++) {
        pos_bits[i] &= other.pos_bits[i];
        neg_bits[i] &= other.neg_bits[i];
    }
    return *this;
}
