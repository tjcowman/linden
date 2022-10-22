
#pragma once

#include "Types.hpp"

namespace Linden::Genetics
{
    struct GenotypeMatrix
    {
        ID_Sample width;
        ID_Snp height;
        std::vector<ID_Genotype> data;

        std::vector<ID_Genotype>::const_iterator rowBegin(size_t row) const
        {
            return data.begin() + row * width;
        }

        std::vector<ID_Genotype>::const_iterator rowEnd(size_t row) const
        {
            return data.begin() + (row+1) * width;
        }
    };
}