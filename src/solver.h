
#pragma once

#include "structures/sat_instance.h"

struct solver {
    inline solver() {};
    virtual void solve() = 0;
    std::string assignment;
};
