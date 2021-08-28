#include <cstdint>
#include <vector>

#include <immintrin.h>

#include "CommonStructs.h"
#include "Types.h"




//Define types and sizes for handling packed genotype representations, can change to facilitate more bitwise parallelism ex 32bit -> 64bit
using PackedGenotype = uint64_t;
constexpr uint8_t PACKED_SIZE = sizeof(PackedGenotype) * 8;
#define POPCOUNT_FUNCTION __builtin_popcountll


ID_Snp popCount(const std::vector<PackedGenotype>& v, ID_Snp begin, ID_Snp distance);

std::vector<PackedGenotype> bitMerge(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2);

ID_Snp popCountAnd(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2, ID_Snp begin, ID_Snp distance);
ID_Snp popCountAnd(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2, ID_Snp begin1, ID_Snp begin2, ID_Snp distance);

ID_Snp popCountAnd_it(std::vector<PackedGenotype>::const_iterator v1, std::vector<PackedGenotype>::const_iterator v2, ID_Snp distance);


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