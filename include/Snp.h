/**
 * @author Tyler Cowman
 * Stores the case and control information for a single SNP, also stores the orginal index in the full dataset.
 * Constructed from plaintext reperesentations of each genotype, i,e. a one byte 0,1,2 per genotype. During
 * construction the genotypes are compressed using an implementation of the BOOST approach from Wan et al.
 */


#ifndef SNP_H
#define SNP_H

#include <iostream>
#include <vector>
#include <array>
#include <stdint.h>
#include <cmath>
#include "CommonStructs.h"
#include "Bitwise.h"

#include "CTable.h"



class Snp
{
    public:
        Snp(ID_Snp index, const GenotypeMatrix& controls, const GenotypeMatrix& cases);
        Snp(const Snp & s1, const Snp & s2); 
        
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
        //float epistasisTest(const Snp & other)const;
        
        //friend std::ostream& operator<< (std::ostream &out, const Snp & snp);

        static void to_serial(std::ostream& os, const Snp& snp);
        static void from_serial(std::istream& is, Snp& snp);
        
        
        static void setDimensions(ID_Sample controls, ID_Sample cases);
    private:
        void packGenotypes( std::vector<uint8_t>::const_iterator begin, std::vector<uint8_t>::const_iterator end, std::vector<uint64_t>& dest);

        //Will consist of controls 0 , 1 , 2 then cases 0 , 1 , 2
        std::vector<PackedGenotype> allSamples_;
       // std::vector<__m256i> samples_;
        ID_Snp index_;
        static SnpDimensions dim;

        //static ID_Sample COW_;
       // static ID_Sample CAW_;
       // static ID_Sample CAS_;
};




#endif //SNP_H
