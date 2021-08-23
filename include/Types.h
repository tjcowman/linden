#pragma once

#include <cstdint>
#include <limits>

using ID_Snp = uint32_t; //Maxium number of input Snps
using ID_Genotype = uint8_t; //Maximum number of genotype values 0, 1, 2 ...
using ID_Sample = uint32_t; //Maximum number of input samples

namespace ID_Invalid{
	constexpr ID_Snp Snp = std::numeric_limits<ID_Snp>::max();
}