#pragma once

#include <cstdint>
#include <limits>
#include <iostream>
#include <vector>

using ID_Snp = uint32_t; //Maxium number of input Snps


using ID_Genotype = uint8_t; //Maximum number of genotype values 0, 1, 2 ...
using ID_Sample = uint32_t; //Maximum number of input samples

namespace ID_Invalid{
	constexpr ID_Snp Snp = std::numeric_limits<ID_Snp>::max();
}




//TODO: put these in some sort of serialization namespace/file
template<class T, class IT>
void vector_to_serial(std::ostream& os, std::vector<T> e) {

	IT length = e.size();
	os.write(reinterpret_cast<const char*>(&length), sizeof(IT));
	//NOTE: The check for length guards against a crash when the there is no data to write (not sure why this crashes)
	if(length >0)
		os.write(reinterpret_cast<const char*>(&e[0]), length * sizeof(T));
}

template<class T, class IT>
std::vector<T> vector_from_serial(std::istream& is) {
	//Read the length
	IT length;
	is.read(reinterpret_cast<char*>(&length), sizeof(IT));

	//Allocate space for the vector
	std::vector<T> e;
	e.resize(length);

	//Read the elements 
	//NOTE: The check for length guards against a crash when the there is no data to write (not sure why this crashes)
	if (length > 0)
		is.read(reinterpret_cast<char*>(&e[0]), length * sizeof(T));

	return e;
}
