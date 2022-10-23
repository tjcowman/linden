
#pragma once

#include <vector>

#include "Locus.hpp"
#include "Snp.hpp"

namespace Linden::Core
{
	class SnpSet
	{
	public:
		// For use with from_serial
		SnpSet() :
			sizeUnfiltered(0),
			dim(Snp::Dimensions(0,0)),
			data(std::vector<Snp>())
		{ }

		// NOTE: CURRENTLY ONLY CONSTRUCT WITH THIS. LOCATIONS MUST BE INDEXED BASED ON THE STORED SNP_INDEX DUE TO FILTERING EX: IN LDFOReet creation
		SnpSet(const std::vector<Genetics::Locus>& loci, const Genetics::GenotypeMatrix& controls, const Genetics::GenotypeMatrix& cases);

		// Getters
		inline Genetics::Id::Snp size() const
		{
			return data.size();
		}

		inline Genetics::Id::Snp getSizeUnfiltered() const
		{
			return sizeUnfiltered;
		}

		inline const Snp::Dimensions& getDimensions() const
		{
			return dim;
		}

		inline const std::vector<Snp>& getSnps() const
		{
			return data;
		}

		inline const std::vector<Genetics::Locus>& getLoci() const
		{
			return loci;
		}

		// Primarily a wrapper for the std library remove_if returns the number of elements removed
		template<typename F>
		Genetics::Id::Snp remove_if(F f)
		{
			auto it = std::remove_if(data.begin(), data.end(), f);
			Genetics::Id::Snp removed = std::distance(it, data.end());
			data.erase(it, data.end());
			return removed;
		}

		inline void truncateTo(Genetics::Id::Snp n)
		{
			data.erase(data.begin()+n, data.end());
		}

		static void to_serial(std::ostream& os, const SnpSet& snpSet);
		static SnpSet from_serial(std::istream& is);


	private:
		std::vector<Genetics::Locus> loci;
		Genetics::Id::Snp sizeUnfiltered;
		Snp::Dimensions dim;
		std::vector<Snp> data;
	};
}