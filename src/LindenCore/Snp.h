/**
 * @author Tyler Cowman
 * Stores the case and control information for a single SNP, also stores the orginal index in the full dataset.
 * Constructed from plaintext reperesentations of each genotype, i,e. a one byte 0,1,2 per genotype. During
 * construction the genotypes are compressed using an implementation of the BOOST approach from Wan et al.
 */


#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <stdint.h>
#include <cmath>
#include <tuple>
#include "CommonStructs.h"
#include "Bitwise.h"

#include "CTable.h"


struct SnpDimensions {

    SnpDimensions() {};

    SnpDimensions(ID_Sample numControls, ID_Sample numCases) :
        numControls_(numControls),
        numCases_(numCases),
        //Determine the number of packed elements required to store num samples, adds one element to handle the last non-full element
        CONR_(numControls / Bitwise::Size + (numControls % Bitwise::Size != 0)),
        CASR_(numCases / Bitwise::Size + (numCases % Bitwise::Size != 0)),
        CASS_(3 * (numControls / Bitwise::Size + (numControls % Bitwise::Size != 0))){}

    ID_Sample numControls_;
    ID_Sample numCases_;
    ID_Sample CONR_;
    ID_Sample CASR_;
    ID_Sample CASS_;
};

class Snp
{
    public:
        Snp();
        Snp(ID_Snp index, const GenotypeMatrix& controls, const GenotypeMatrix& cases);
        Snp(const Snp & s1, const Snp & s2); 

        bool operator==(const Snp& lhs)const;
        
        //Getters
        ID_Snp getIndex()const;
        
        //Calculations
        float computeMinorAlleleFrequency()const;
        
        //Comparisons
        ID_Sample computeDifferences(const Snp & other)const;
        float computeUnknownRatio()const;
        
        //Epistasis testing
        static void fillTable(CTable2& t, const Snp& snp1, const Snp& snp2);

        float marginalTest()const;


        static void to_serial(std::ostream& os, const Snp& snp);
        static Snp from_serial(std::istream& is);
        
        
        static void setDimensions(ID_Sample controls, ID_Sample cases);
        static const SnpDimensions& getDimensions(); //Mostly used as access for serialization/deserialization

    private:
        void packGenotypes( std::vector<uint8_t>::const_iterator begin, std::vector<uint8_t>::const_iterator end, std::vector<uint64_t>& dest);

        //Will consist of controls 0 , 1 , 2 then cases 0 , 1 , 2
        std::vector<Bitwise::Genotype> allSamples_;
       // std::vector<__m256i> samples_;
        ID_Snp index_;
        static SnpDimensions dim;

        //static ID_Sample COW_;
       // static ID_Sample CAW_;
       // static ID_Sample CAS_;
};
