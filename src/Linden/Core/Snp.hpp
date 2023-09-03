
#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <stdint.h>
#include <tuple>
#include <vector>

#include "Bitwise.h"
#include "ContingencyTable.hpp"
#include "GenotypeMatrix.hpp"

namespace Linden::Core
{
    ////////////////////////////////////////////////////////////////////////////////
    //! @class Snp
    //!
    //! Stores the case and control information for a single SNP, also stores the 
    //! orginal index in the full dataset. Constructed from plaintext 
    //! reperesentations of each genotype, i,e. a one byte 0,1,2 per genotype.
    //! During construction the genotypes are compressed using an implementation of
    //! the BOOST approach from Wan et al.
    ////////////////////////////////////////////////////////////////////////////////
    class Snp
    {
    public:
        ////////////////////////////////////////////////////////////////////////////
        //! @struct Dimensions
        //!
        //! Stores data relating to the size input snp data including conversions
        //! to the number of packed elements and their locations in the underlying 
        //! packed genotype array.
        ////////////////////////////////////////////////////////////////////////////
        struct Dimensions
        {
            Dimensions()
            { }

            Dimensions(Genetics::Id::Sample numControls, Genetics::Id::Sample numCases) :
                numControls_(numControls),
                numCases_(numCases),
                // Determine the number of packed elements required to store num samples
                // adds one element to handle the last non-full element
                numPackedControls_(numControls / Bitwise::Size + (numControls % Bitwise::Size != 0)),
                numPackedCases_(numCases / Bitwise::Size + (numCases % Bitwise::Size != 0)),
                casesBegin_(3 * (numControls / Bitwise::Size + (numControls % Bitwise::Size != 0)))
            { }

            //! Total number of control samples.
            Genetics::Id::Sample numControls_;
            //! Total number of case samples.
            Genetics::Id::Sample numCases_;
            //! Number of packed control elements.
            Genetics::Id::Sample numPackedControls_;
            //! Number of packed case elements.
            Genetics::Id::Sample numPackedCases_;
            //! Starting index of the packed case samples.
            Genetics::Id::Sample casesBegin_;
        };

        // Construct a Snp from an index and vector of samples, ex: from serialized data
        Snp(Genetics::Id::Snp index, const std::vector<Bitwise::Genotype>& samples);

        // Construct a Snp from an index in a case and control matrix.
        Snp(Genetics::Id::Snp index, const Genetics::GenotypeMatrix& controls, const Genetics::GenotypeMatrix& cases);
        
        // Constructs a new snp by merging two children snps. Note that the
        Snp(const Snp & s1, const Snp & s2);

        // Equivalence operator compares index and genotype values.
        inline bool operator==(const Snp& lhs) const
        {
            return std::tie(index_, allSamples_) == std::tie(lhs.index_, lhs.allSamples_);
        }

        // Gets the unique index of the snp.
        inline Genetics::Id::Snp getIndex() const
        {
            return index_;
        }

        // Sets the snp dimensions for all snps.
        inline static void setDimensions(Genetics::Id::Sample controls, Genetics::Id::Sample cases)
        {
            Snp::dim = Dimensions(controls, cases);
        }
        // Gets the snp dimensions.
        inline static const Dimensions& getDimensions()
        {
            return Snp::dim;
        };

        // Computes the frequency of the less common allele.
        float computeMinorAlleleFrequency() const;

        // Calculate the number of differnent genotypes from another snp.
        Genetics::Id::Sample computeDifferences(const Snp & other) const;

        // Calculate proportion of ambiguous genotypes.
        float computeUnknownRatio() const;

        // Fills out a contingency table with genotypes from two snps.
        static void fillTable(Statistics::ContingencyTable<9>& t, const Snp& snp1, const Snp& snp2);

        // Peforms a marginal significance test on the snp.
        float marginalTest() const;

        static void to_serial(std::ostream& os, const Snp& snp);
        static Snp from_serial(std::istream& is);

    private:

        // Gets iterator to the start of packed control genotpyes
        inline std::vector<Bitwise::Genotype>::const_iterator GetControlsBegin(int genotype) const
        {
            return allSamples_.begin() + (genotype * dim.numPackedControls_);
        }
        // Gets iterator to the start of packed case genotpyes
        inline std::vector<Bitwise::Genotype>::const_iterator GetCasesBegin(int genotype) const
        {
            return allSamples_.begin() + dim.casesBegin_ + (genotype * dim.numPackedCases_);
        }

        // Packs range of genotype values into their bitwise representations.
        void packGenotypes(std::vector<uint8_t>::const_iterator begin, 
                        std::vector<uint8_t>::const_iterator end, 
                        std::vector<Bitwise::Genotype>& dest);

        //! Unique index corresponding to the snp.
        Genetics::Id::Snp index_;
        //! Packed Genotypes controls 0,1,2 then cases 0,1,2.
        std::vector<Bitwise::Genotype> allSamples_;

        //! Dimensions of the current snp data
        static Dimensions dim;
    };
}