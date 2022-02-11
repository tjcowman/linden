#pragma once


#include "Snp.h"
#include "CommonStructs.h"

#include <vector>


class SnpSet {
public:

	SnpSet(); //For use with from_serial
	//NOTE: CURRENTLY ONLY CONSTRUCT WITH THIS. LOCATIONS MUST BE INDEXED BASED ON THE STORED SNP_INDEX DUE TO FILTERING EX: IN LDFOReet creation
	SnpSet(const std::vector<Locus>& loci, const GenotypeMatrix& controls, const GenotypeMatrix& cases);

	ID_Snp size()const;
	ID_Snp getSizeUnfiltered()const;
	const SnpDimensions& getDimensions()const;

	const std::vector<Snp>& getSnps();
	//const std::vector<Location>& getLocations();

	//Primarily a wrapper for the std library remove_if returns the number of elements removed
	template<typename F>
	ID_Snp remove_if(F f);
	
	void truncateTo(ID_Snp n)
	{
		data.resize(n);
	}


	static void to_serial(std::ostream& os, const SnpSet& snpSet);
	static SnpSet from_serial(std::istream& is);


	std::vector<Locus> loci;
private:
	ID_Snp sizeUnfiltered;
	SnpDimensions dim;
	std::vector<Snp> data;
};

template<typename F>
ID_Snp SnpSet::remove_if(F f) {

	auto it = std::remove_if(data.begin(), data.end(), f);
	ID_Snp removed = std::distance(it, data.end());
	data.erase(it, data.end());
	return removed;
}