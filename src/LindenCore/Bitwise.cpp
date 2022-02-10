#include "Bitwise.h"


ID_Snp popCount(const std::vector<PackedGenotype>& v, ID_Snp begin, ID_Snp distance){
    PackedGenotype retVal = 0;
    for (ID_Snp i = begin; i < begin + distance; ++i)
        retVal += static_cast<PackedGenotype>(POPCOUNT_FUNCTION(v[i]));

    return static_cast<ID_Snp>(retVal);
}

std::vector<PackedGenotype> bitMerge(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2){
    std::vector<PackedGenotype>  retVal;
    retVal.reserve(v1.size());

    for (ID_Snp i = 0; i < v1.size(); ++i)
        retVal.push_back(v1[i] & v2[i]);

    return retVal;
}

ID_Snp popCountAnd(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2, ID_Snp begin, ID_Snp distance){
    PackedGenotype retVal = 0;

    for (ID_Snp i = begin; i < begin + distance; ++i)
    {
        if ((v1[begin] | v2[begin]) != 0)
            retVal += static_cast<PackedGenotype>(POPCOUNT_FUNCTION(v1[i] & v2[i]));
    }
    return static_cast<ID_Snp>(retVal);
}

ID_Snp popCountAnd(const std::vector<PackedGenotype>& v1, const std::vector<PackedGenotype>& v2, ID_Snp begin1, ID_Snp begin2, ID_Snp distance){
    PackedGenotype retVal = 0;
    for (ID_Snp i = 0; i < distance; ++i)
    {
        //There are many instances of long strings of zeroes, primarily due to homozygous minor this saves popcount operations
        if ((v1[begin1 + i] | v2[begin2 + i]) != 0)
            retVal += static_cast<PackedGenotype>(POPCOUNT_FUNCTION(v1[begin1 + i] & v2[begin2 + i]));
    }
    return static_cast<ID_Snp>(retVal);
}

ID_Snp popCountAnd_it(std::vector<PackedGenotype>::const_iterator v1, std::vector<PackedGenotype>::const_iterator v2, ID_Snp distance) {
    PackedGenotype retVal = 0;

    for (size_t i = 0; i < distance; ++i) {
        if ((*v1 | *v2) != 0) {
            retVal += static_cast<PackedGenotype>(POPCOUNT_FUNCTION(*v1 & *v2));
        }
        ++v1;
        ++v2;
    }
    return static_cast<ID_Snp>(retVal);
}