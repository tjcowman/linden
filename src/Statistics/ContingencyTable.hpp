
#pragma once

#include <array>
#include <cmath>

#include "CommonStructs.h"

namespace Statistics
{
    //Used for determining the chi2 correseponding to a pvalue of 1*10^-i for single locus significance
    //Input as an integer denoting the -log10 signficance up to 6
    //inline const static std::array<float, 7> chi2DegreesFreedomTable{0.0f, 2.71f, 6.63f, 10.82f, 15.14f, 19.57f, 24.87f};
    constexpr auto CreateChi2Table()
    {
        auto table = std::array<float, 7>{0.0f, 2.71f, 6.63f, 10.82f, 15.14f, 19.57f, 24.87f};
        return table;
    }
    constexpr auto chi2Table = CreateChi2Table();
    inline const std::array<float, 7>& GetChi2Table()
    {
        return chi2Table;
    }
}

template<std::size_t Values>
class ContingencyTable
{
public:

    inline void zero()
    {
        m_data.fill(0);
    }

    inline std::array<ID_Sample, Values * 3 + 2>& Data()
    {
        return m_data;
    }

    ID_Sample& ColTotal(int col)
    {
        return m_data[Values*3 + col];
    }

    ID_Sample& RowTotal(int row)
    {
        return m_data[(row*3) + 2];
    }

    inline float Chi2()
    {
        float chi2 = 0.0;
        
        auto c3 = m_data[colTotalsIndex] + m_data[colTotalsIndex+1]; // total
        for (int row = 0; row < Values; ++row)
        {
            for (int col = 0; col < levels; ++col)
            {
                ID_Sample c1 = m_data[row * 3 + levels]; // row total
                ID_Sample c2 = m_data[colTotalsIndex + col]; // col total
                ID_Sample c4 = m_data[row * 3 + col];

                float expected = (float)(c1 * c2) / (float)c3;
                chi2 += (pow(c4 - expected, 2) / expected);
            }
        }
        return chi2;
    }

private:
    //! Starting index of the column totals.
    const static int colTotalsIndex = Values * 3;
    //! Number of levels ex: case, control groups
    const static int levels = 2;
    //! Backing data structure.
    std::array<ID_Sample, Values * 3 + 2> m_data;
};

namespace ContingencyTableTmp
{
    const static std::array<std::pair<uint8_t, uint8_t>, 9> rowOrder
    {{ 
        {0,0}, {0,1}, {0,2},
        {1,0}, {1,1}, {1,2},
        {2,0}, {2,1}, {2,2}
    }};
}
