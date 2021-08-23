#pragma once

#include "Types.h"
#include <array>

#include "Snp.h"
class Snp;

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

    float chi2() {
        float chi2 = 0.0;
        constexpr uint8_t factors = 3;
        constexpr uint8_t levels = 2;

        for (uint8_t row = 0; row < factors * factors; ++row) {
            for (uint8_t col = 0; col < levels; ++col) {
                ID_Sample c1 = data[row * 3 + 2];
                ID_Sample c2 = data[27 + col];
                ID_Sample c3 = data[27] + data[28];
                ID_Sample c4 = data[row * 3 + col];

                float expected = (float)(c1 * c2) / (float)c3;
                chi2 += (pow(c4 - expected, 2) / expected);
            }
        }

        return chi2;
    }


};