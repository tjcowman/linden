
#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <stdint.h>
#include <cmath>
#include <tuple>
#include "Bitwise.h"

#include "ContingencyTable.hpp"

struct SnpDimensions
{
    SnpDimensions()
    { }

    SnpDimensions(ID_Sample numControls, ID_Sample numCases) :
        numControls_(numControls),
        numCases_(numCases),
        //Determine the number of packed elements required to store num samples, adds one element to handle the last non-full element
        CONR_(numControls / Bitwise::Size + (numControls % Bitwise::Size != 0)),
        CASR_(numCases / Bitwise::Size + (numCases % Bitwise::Size != 0)),
        CASS_(3 * (numControls / Bitwise::Size + (numControls % Bitwise::Size != 0)))
    { }

    ID_Sample numControls_;
    ID_Sample numCases_;
    ID_Sample CONR_;
    ID_Sample CASR_;
    ID_Sample CASS_;
};

////////////////////////////////////////////////////////////////////////////////
//! @class Snp
//! Stores the case and control information for a single SNP, also stores the 
//! orginal index in the full dataset. Constructed from plaintext 
//! reperesentations of each genotype, i,e. a one byte 0,1,2 per genotype.
//! During construction the genotypes are compressed using an implementation of
//! the BOOST approach from Wan et al.
////////////////////////////////////////////////////////////////////////////////
class Snp
{
public:
    // Construct a Snp from an index and vector of samples, ex: from serialized data
    Snp(ID_Snp index, const std::vector<Bitwise::Genotype>& samples);

    // Construct a Snp from an index in a case and control matrix.
    Snp(ID_Snp index, const GenotypeMatrix& controls, const GenotypeMatrix& cases);
    
    // Constructs a new snp by merging two children snps. Note that the
    Snp(const Snp & s1, const Snp & s2);

    // Equivalence operator compares index and genotype values.
    inline bool operator==(const Snp& lhs) const
    {
        return std::tie(index_, allSamples_) == std::tie(lhs.index_, lhs.allSamples_);
    }

    // Gets the unique index of the snp.
    inline ID_Snp getIndex() const
    {
        return index_;
    }

    // Sets the snp dimensions for all snps.
    inline static void setDimensions(ID_Sample controls, ID_Sample cases)
    {
        Snp::dim = SnpDimensions(controls, cases);
    }
    // Gets the snp dimensions.
    inline static const SnpDimensions& getDimensions()
    {
        return Snp::dim;
    };

    // Computes the frequency of the less common allele.
    float computeMinorAlleleFrequency() const;

    // Calculate the number of differnent genotypes from another snp.
    ID_Sample computeDifferences(const Snp & other) const;

    // Calculate proportion of ambiguous genotypes.
    float computeUnknownRatio()const;

    // Fills out a contingency table with genotypes from two snps.
    static void fillTable(ContingencyTable2& t, const Snp& snp1, const Snp& snp2);

    // Peforms a marginal significance test on the snp.
    float marginalTest() const;

    static void to_serial(std::ostream& os, const Snp& snp);
    static Snp from_serial(std::istream& is);

private:
    // Packs range of genotype values into their bitwise representations.
    void packGenotypes(std::vector<uint8_t>::const_iterator begin, 
        std::vector<uint8_t>::const_iterator end, 
        std::vector<Bitwise::Genotype>& dest);

    //! Unique index corresponding to the snp.
    ID_Snp index_;
    //! Packed Genotypes controls 0,1,2 then cases 0,1,2.
    std::vector<Bitwise::Genotype> allSamples_;

    //! Dimensions of the current snp data
    static SnpDimensions dim;
};
