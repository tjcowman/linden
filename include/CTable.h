#pragma once

#define GENOTYPE_LEVELS 3
#define GENOTYPE_PAIRINGS 9
#define CONTINGENCY_COLUMNS 2

struct CTable1 {
    std::array<int, 6> M_;
    std::array<int, 3> RT_;
    std::array<int, 2> CT_;

    int& a(int x, int y)
    {
        return M_[x + (3 * y)];
    }

    void addOne()
    {
        for (int i = 0; i < M_.size(); ++i)
            M_[i]++;
    }

    void calculateTotals()
    {
        CT_.fill(0);
        RT_.fill(0);

        for (int i = 0; i < 3; ++i)
        {
            RT_[i] = M_[i] + M_[i + 3];
            CT_[0] += M_[i];
            CT_[1] += M_[i + 3];
        }
    }
};

struct CTable2 {
    std::array<int, 29> M_;

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