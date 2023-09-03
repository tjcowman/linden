
#pragma once

#include <limits>

namespace Genetics::Id
{
    // Maxium number of input Snps
    using Snp = std::uint32_t;
    // Maximum number of genotype values 0, 1, 2 ...
    using Genotype = std::uint8_t;
    // Maximum number of input samples
    using Sample = std::uint32_t; 

    constexpr Snp InvalidSnp = std::numeric_limits<Snp>::max();
}