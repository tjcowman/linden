
#pragma once

#include "Id.hpp"

namespace Genetics
{
    struct GenotypeMatrix
    {
        Genetics::Id::Sample width;
        Genetics::Id::Snp height;
        std::vector<Genetics::Id::Genotype> data;

        std::vector<Genetics::Id::Genotype>::const_iterator rowBegin(size_t row) const
        {
            return data.begin() + row * width;
        }

        std::vector<Genetics::Id::Genotype>::const_iterator rowEnd(size_t row) const
        {
            return data.begin() + (row+1) * width;
        }
    };
}