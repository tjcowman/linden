#pragma once


#include "Snp.h"
#include "CommonStructs.h"
#include "Types.h"
#include <vector>


class SnpSet {
public:
	SnpSet(const std::vector<Locus>& loci, const GenotypeMatrix& controls, const GenotypeMatrix& cases);

	const std::vector<Snp>& getSnps();
	const std::vector<Location>& getLocations();

	//Primarily a wrapper for the std library remove_if returns the number of elements removed
	template<typename F>
	ID_Snp remove_if(F f);
	
	//TODO SET SNP DIM

private:
	ID_Snp sizeUnfiltered;
	SnpDimensions dim;
	std::vector<Snp> data;
	std::vector<Location> locations;
};

template<typename F>
ID_Snp SnpSet::remove_if(F f) {

	auto it = std::remove_if(data.begin(), data.end(), f);
	ID_Snp removed = std::distance(it, data.end());
	data.erase(it, data.end());
	return removed;
}