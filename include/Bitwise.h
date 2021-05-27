#include <cstdint>
#include <vector>

#include <immintrin.h>

#include "CommonStructs.h"
#include "Types.h"




//Define types and sizes for handling packed genotype representations, can change to facilitate more bitwise parallelism ex 32bit -> 64bit
using PackedGenotype = uint64_t;
constexpr uint8_t PACKED_SIZE = sizeof(PackedGenotype) * 8;

#define POPCOUNT_FUNCTION __builtin_popcountll

static ID_Snp popCount(const std::vector<PackedGenotype>& v, ID_Snp begin, ID_Snp distance)
{
    PackedGenotype retVal = 0;
    for (ID_Snp i = begin; i < begin + distance; ++i)
        retVal += static_cast<PackedGenotype>(POPCOUNT_FUNCTION(v[i]));

    return static_cast<ID_Snp>(retVal);
}

static std::vector<PackedGenotype> bitMerge(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2)
{
    std::vector<PackedGenotype>  retVal;
    retVal.reserve(v1.size());

    for (ID_Snp i = 0; i < v1.size(); ++i)
        retVal.push_back(v1[i] & v2[i]);

    return retVal;
}

static ID_Snp popCountAnd(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2, ID_Snp begin, ID_Snp distance)
{
    PackedGenotype retVal = 0;

    for (ID_Snp i = begin; i < begin + distance; ++i)
    {
        if ((v1[begin] | v2[begin]) != 0)
            retVal += static_cast<PackedGenotype>(POPCOUNT_FUNCTION(v1[i] & v2[i]));
    }
    return static_cast<ID_Snp>(retVal);
}

static ID_Snp popCountAnd(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2, ID_Snp begin1, ID_Snp begin2, ID_Snp distance)
{
    PackedGenotype retVal = 0;
    for (ID_Snp i = 0; i < distance; ++i)
    {
        //There are many instances of long strings of zeroes, primarily due to homozygous minor this saves popcount operations
        if ((v1[begin1 + i] | v2[begin2 + i]) != 0)
            retVal += static_cast<PackedGenotype>(POPCOUNT_FUNCTION(v1[begin1 + i] & v2[begin2 + i]));
    }
    return static_cast<ID_Snp>(retVal);
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