#pragma once

#include "Types.h"

#include<array>

#define GENOTYPE_LEVELS 3
#define GENOTYPE_PAIRINGS 9
#define CONTINGENCY_COLUMNS 2

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

    int& a(int x, int y) {
        return M_[x + (9 * y)];
    }

    void printM() {
        std::cout << M_[0] << " " << M_[1] << "\n"
            << M_[2] << " " << M_[3] << "\n"
            << M_[4] << " " << M_[5] << "\n"
            << M_[6] << " " << M_[7] << "\n"
            << M_[8] << " " << M_[9] << "\n"
            << M_[10] << " " << M_[11] << "\n"
            << M_[12] << " " << M_[13] << "\n"
            << M_[14] << " " << M_[15] << "\n"
            << M_[16] << " " << M_[17] << "\n"
            << "\n";
    }


    void calculateTotals() {
        for (int i = 0; i < 9; ++i) {
            M_[18 + i] = M_[i] + M_[i + 9];
            M_[27] += M_[i];
            M_[28] += M_[i + 9];
        }
    }
};