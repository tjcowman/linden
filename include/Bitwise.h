#include <cstdint>
#include <vector>
#include "CommonStructs.h"
#include "Types.h"


//Define types and sizes for handling packed genotype representations, can change to facilitate more bitwise parallelism ex 32bit -> 64bit
using PackedGenotype = uint64_t;
constexpr uint8_t PACKED_SIZE = sizeof(PackedGenotype) * 8;


#include <intrin.h> 

#define POPCOUNT_FUNCTION __builtin_popcountll



static PackedGenotype popCount(const std::vector<PackedGenotype>& v, ID_Snp begin, ID_Snp distance)
{
    PackedGenotype retVal = 0;
    for (ID_Snp i = begin; i < begin + distance; ++i)
        retVal += POPCOUNT_FUNCTION(v[i]);

    return retVal;
}

static std::vector<PackedGenotype> bitMerge(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2)
{
    std::vector<PackedGenotype>  retVal;
    retVal.reserve(v1.size());

    for (ID_Snp i = 0; i < v1.size(); ++i)
        retVal.push_back(v1[i] & v2[i]);

    return retVal;
}

static PackedGenotype popCountAnd(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2, ID_Snp begin, ID_Snp distance)
{
    PackedGenotype retVal = 0;

    for (ID_Snp i = begin; i < begin + distance; ++i)
    {
        if ((v1[begin] | v2[begin]) != 0)
            retVal += POPCOUNT_FUNCTION(v1[i] & v2[i]);
    }
    return retVal;
}

static PackedGenotype popCountAnd(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2, ID_Snp begin1, ID_Snp begin2, ID_Snp distance)
{
    PackedGenotype retVal = 0;
    //#pragma acc parallel reduction(+:retVal)
    for (ID_Snp i = 0; i < distance; ++i)
    {
        //There are many instnaces of long strings of zeroes, primarily due to homozygous minor this saves popcount operations
        if ((v1[begin1 + i] | v2[begin2 + i]) != 0)
            retVal += POPCOUNT_FUNCTION(v1[begin1 + i] & v2[begin2 + i]);
    }
    return retVal;
}