
#pragma once

#include "Snp.hpp"
#include "CommonStructs.h"

#include <vector>

class SnpSet
{
public:
	 // For use with from_serial
	SnpSet() :
		sizeUnfiltered(0),
		dim(SnpDimensions(0,0)),
		data(std::vector<Snp>())
	{ }

	// NOTE: CURRENTLY ONLY CONSTRUCT WITH THIS. LOCATIONS MUST BE INDEXED BASED ON THE STORED SNP_INDEX DUE TO FILTERING EX: IN LDFOReet creation
	SnpSet(const std::vector<Locus>& loci, const GenotypeMatrix& controls, const GenotypeMatrix& cases);

	// Getters
	inline ID_Snp size() const
	{
		return data.size();
	}

	inline ID_Snp getSizeUnfiltered() const
	{
		return sizeUnfiltered;
	}

	inline const SnpDimensions& getDimensions() const
	{
		return dim;
	}

	inline const std::vector<Snp>& getSnps() const
	{
		return data;
	}

	inline const std::vector<Locus>& getLoci() const
	{
		return loci;
	}

	// Primarily a wrapper for the std library remove_if returns the number of elements removed
	template<typename F>
	ID_Snp remove_if(F f)
	{
		auto it = std::remove_if(data.begin(), data.end(), f);
		ID_Snp removed = std::distance(it, data.end());
		data.erase(it, data.end());
		return removed;
	}

	inline void truncateTo(ID_Snp n)
	{
		data.erase(data.begin()+n, data.end());
	}

	static void to_serial(std::ostream& os, const SnpSet& snpSet);
	static SnpSet from_serial(std::istream& is);


private:
	std::vector<Locus> loci;
	ID_Snp sizeUnfiltered;
	SnpDimensions dim;
	std::vector<Snp> data;
};
