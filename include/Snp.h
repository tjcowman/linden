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


#define GENOTYPE_LEVELS 3
#define GENOTYPE_PAIRINGS 9
#define CONTINGENCY_COLUMNS 2

#define PACK_SIZE 64
#define PACK_TYPE uint64_t

//If windows use it's popcount



#ifdef _WIN32
    #include <intrin.h> 
    #define POPCOUNT_FUNCTION __popcnt64
#else
    #define POPCOUNT_FUNCTION __builtin_popcountll
#endif





static int popCount(const std::vector<PACK_TYPE> & v, int begin, int distance)
{
    int retVal = 0;
    for(int i=begin; i<begin+distance; ++i)
        retVal += POPCOUNT_FUNCTION(v[i]);

    return retVal;
}

static std::vector<PACK_TYPE> bitMerge(const std::vector<PACK_TYPE> & v1, const std::vector<PACK_TYPE> & v2)
{
    std::vector<PACK_TYPE>  retVal;
    retVal.reserve(v1.size());
    
    for(int i=0; i<v1.size(); ++i)
        retVal.push_back(v1[i] & v2[i]);
    
    return retVal;
}

static int popCountAnd(const std::vector<PACK_TYPE> & v1, const std::vector<PACK_TYPE> & v2, int begin, int distance)
{
    int retVal = 0;

    for(int i=begin; i<begin+distance; ++i)
    {
        if(v1[begin] | v2[begin] !=0)
        retVal += POPCOUNT_FUNCTION(v1[i] & v2[i]);
    }
    return retVal;
}

static int popCountAnd(const std::vector<PACK_TYPE> & v1, const std::vector<PACK_TYPE> & v2, int begin1, int begin2, int distance)
{
    int retVal = 0;
    //#pragma acc parallel reduction(+:retVal)
    for(int i=0; i<distance; ++i)
    {
        //There are many instnaces of long strings of zeroes, primarily due to homozygous minor this saves popcount operations
        if(v1[begin1+i] | v2[begin2+i] !=0)
            retVal += POPCOUNT_FUNCTION(v1[begin1 +i] & v2[begin2 +i]);
    }
    return retVal;
}

struct SmallContingencyTable
{
    std::array<int,6> M_;
    std::array<int,3> RT_;
    std::array<int,2> CT_;
    
    int & a(int x, int y)
    {
        return M_[x+(3*y)];
    }
    
    void addOne()
    {
        for(int i=0; i<M_.size(); ++i)
            M_[i]++;
    }
    
    void calculateTotals()
    {
        CT_.fill(0);
        RT_.fill(0);
        
        for(int i=0; i<3; ++i)
        {
            RT_[i] = M_[i] + M_[i+3];
            CT_[0] += M_[i];
            CT_[1] += M_[i+3];
        }
    }
    
};

struct ContingencyTable
{
    std::array<int,29> M_;

    int & a(int x, int y)
    {
        return M_[x+(9*y)];
    }
    
    void printM()
    {
            
        std::cout<<M_[0]<<" "<<M_[1]<<"\n"
        <<M_[2]<<" "<<M_[3]<<"\n"
        <<M_[4]<<" "<<M_[5]<<"\n"
        <<M_[6]<<" "<<M_[7]<<"\n"
        <<M_[8]<<" "<<M_[9]<<"\n"
        <<M_[10]<<" "<<M_[11]<<"\n"
        <<M_[12]<<" "<<M_[13]<<"\n"
        <<M_[14]<<" "<<M_[15]<<"\n"
        <<M_[16]<<" "<<M_[17]<<"\n"
        <<"\n";
    }
    
    
    void calculateTotals()
    {

        for(int i=0; i<9; ++i)
        {

            M_[18+i] = M_[i] + M_[i+9];
            M_[27] += M_[i];
            M_[28] += M_[i+9];
        }
        
    }
};



class Snp
{
    public:
        Snp(int index, const std::vector<char> & controls, const std::vector<char> & cases);
        Snp(const Snp & cpy);
        
        Snp(const Snp & s1, const Snp & s2); 
        
        //Getters
        int getIndex()const;
        
        //Calculations
        float computeMinorAlleleFrequency()const;
        
        //Comparisons
        int computeDifferences(const Snp & other)const;
        float computeUnknownRatio()const;
        
        //Epistasis testing
        float marginalTest()const;
        float epistasisTest(const Snp & other)const;
        
        friend std::ostream& operator<< (std::ostream &out, const Snp & snp);
        
    private:
        
        std::array<std::vector<PACK_TYPE>, GENOTYPE_LEVELS> packGenotypesHelper(const std::vector<char> & genotypes); 
        void packGenotypes(std::vector<PACK_TYPE> & allPackedGenotypes, const std::vector<char> & controlGenotypes, const std::vector<char> & caseGenotypes);

        //Will consist of controls 0 , 1 , 2 then cases 0 , 1 , 2
        std::vector<PACK_TYPE> allSamples_;

        int index_;
        
        static int numCases_;
        static int numControls_;
        static int CONR_;
        static int CASR_;
        static int CASS_;
        
        
        
};




#endif //SNP_H
