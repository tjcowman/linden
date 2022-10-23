
#pragma once

#include <cstdint>
#include <vector>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

#include "Id.hpp"

//Define types and sizes for handling packed genotype representations, can change to facilitate more bitwise parallelism ex 32bit -> 64bit
namespace Linden::Bitwise
{
    // Compiler specific function for performing population counting
#if defined(_MSC_VER)
    #define POPCOUNT_FUNCTION __popcnt
        // Underlying element type to store consecutive genotype bits in
    using Genotype = std::uint32_t;
#else
    #define POPCOUNT_FUNCTION __builtin_popcountll
    // Underlying element type to store consecutive genotype bits in
    using Genotype = std::uint64_t;
#endif

    // Number of bits in one primitive element
    constexpr uint8_t Size = sizeof(Genotype) * 8;

    // Produces a new vector of Genotypes by peforming a bitwise and
    inline std::vector<Genotype> merge(const std::vector<Genotype>& v1, const std::vector<Genotype>& v2)
    {
        std::vector<Genotype>  retVal;
        retVal.reserve(v1.size());

        for (Genetics::Id::Snp i = 0; i < v1.size(); ++i)
        {
            retVal.push_back(v1[i] & v2[i]);
        }
        return retVal;
    }

    // Peforms population count on a vector of genotype data
    inline Genetics::Id::Snp count(std::vector<Genotype>::const_iterator  v, Genetics::Id::Snp distance)
    {
        Genotype retVal = 0;
        for (Genetics::Id::Snp i = 0; i < distance; ++i)
        {
            retVal += static_cast<Genotype>(POPCOUNT_FUNCTION(*v));
            ++v;
        }
        return static_cast<Genetics::Id::Snp>(retVal);
    }

    // Peforms bitwise and followdby population count on two parallel vectors of genotype data
    inline Genetics::Id::Snp andCount(std::vector<Genotype>::const_iterator v1, std::vector<Genotype>::const_iterator v2, Genetics::Id::Snp distance)
    {
        Genotype retVal = 0;

        for (size_t i = 0; i < distance; ++i)
        {
            retVal += static_cast<Genotype>(POPCOUNT_FUNCTION(*v1 & *v2));

            ++v1;
            ++v2;
        }
        return static_cast<Genetics::Id::Snp>(retVal);
    }


    //SIMD TEST
    /*
    static PackedGenotype popCountAnd_Vec(const std::vector<__m256i>& v1, const std::vector<__m256i>& v2, ID_Snp begin1, ID_Snp begin2, ID_Snp distance)
    {
        PackedGenotype retVal = 0;
        __m256i res;
        for (ID_Snp i = 0; i < distance; ++i)
        {
            res = _mm256_and_si256(v1[begin1 +i], v2[begin2 +i]);

            if(res[0] != 0)
            retVal += POPCOUNT_FUNCTION(res[0]);
            if (res[1] != 0)
            retVal += POPCOUNT_FUNCTION(res[1]);
            if (res[2] != 0)
            retVal += POPCOUNT_FUNCTION(res[2]);
            if (res[3] != 0)
            retVal += POPCOUNT_FUNCTION(res[3]);

        }
        //std::cout << "-------------------------\n";
        return retVal;
    }
    */
}