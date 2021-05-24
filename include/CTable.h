#pragma once

#include "Types.h"

#include<array>

//#define GENOTYPE_LEVELS 3
//#define GENOTYPE_PAIRINGS 9
//#define CONTINGENCY_COLUMNS 2

namespace CTable {
    const static std::array<std::pair<uint8_t, uint8_t>, 9> rowOrder{ { {0,0},{0,1},{0,2},{1,0},{1,1},{1,2},{2,0},{2,1},{2,2} } };
}


struct CTable1 {
    std::array<ID_Sample, 11> data;

    void zero(){
        data.fill(0);
    }
};

struct CTable2 {
    std::array<int, 29> M_;
    std::array<ID_Sample, 29> data;

    void zero() {
        data.fill(0);
    }
};