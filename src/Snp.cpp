#include "Snp.h"

//Define static members, will be set correctly in the Snp constructors
int Snp::numControls_;
int Snp::numCases_;
int Snp::CONR_;
int Snp::CASR_;
int Snp::CASS_;


Snp::Snp(int index, const std::vector<char> & controls, const std::vector<char> & cases)
{

    packGenotypes(allSamples_, controls, cases);
    
    index_ = index;
    
    numControls_ = controls.size();
    numCases_ = cases.size();
    CONR_ = ceil(numControls_/(double)PACK_SIZE);
    CASR_ = ceil(numCases_/(double)PACK_SIZE);
    CASS_ = (3*CONR_) ;//+ 1;
    
}

Snp::Snp(const Snp & cpy)
{
    index_ = cpy.index_;
    allSamples_ = cpy.allSamples_;
}

Snp::Snp(const Snp & s1, const Snp & s2)
{
    //Index is -1 because any node created this way will be internal
    //and not refer to a specific snp from the dataset
    index_ = -1;
    allSamples_ = bitMerge(s1.allSamples_, s2.allSamples_);
}

int Snp::getIndex()const
{
    return index_;
}


float Snp::computeMinorAlleleFrequency()const
{
    int homoMajor = popCount(allSamples_, 0, CONR_) + popCount(allSamples_, CASS_, CASR_) ;
    int hetero = popCount(allSamples_, CONR_, CONR_) + popCount(allSamples_, CASS_+CASR_, CASR_); 
    int homoMinor = popCount(allSamples_, 2*CONR_, CONR_) + popCount(allSamples_, CASS_+(2*CASR_), CASR_) ;
    
    return (2*homoMinor+hetero)/(float)(2*homoMajor + hetero + 2*homoMinor);
}


int Snp::computeDifferences(const Snp & other)const
{
    return ((numControls_ + numCases_)- popCountAnd(allSamples_, other.allSamples_, 0, allSamples_.size()));
}


float Snp::computeUnknownRatio()const
{
    int numberKnown = popCount(allSamples_, 0, allSamples_.size());

    return ( numberKnown/(float)(numControls_ + numCases_) );
}

float Snp::marginalTest()const
{
    float retVal=0.0;
    SmallContingencyTable t;
    
    for(int i=0; i<GENOTYPE_LEVELS; ++i)
    {
        t.M_[i] = popCount(allSamples_, i*CONR_, CONR_);
        t.M_[i+GENOTYPE_LEVELS] = popCount(allSamples_, CASS_ + (i*CASR_), CASR_);
    }
       
    t.addOne(); //For correction    
    t.calculateTotals();
        
		
    int factors = 3;
    int levels = 2;
    
    for(int y=0; y< factors; ++y)
    {

        for(int x=0; x< levels; ++x)
        {            
            //Use margin values
            float c1 = (float)t.RT_[y];
            float c2 = (float)t.CT_[x];
            float c3 = (float)(t.CT_[0] + t.CT_[1]);          
            
            float c4 = (float)t.a(y,x);

            float expected  = ( c1 * c2 ) / c3 ;

            retVal = retVal +  (pow(c4-expected,2)/expected);

        }
    }
    return retVal;
}

float Snp::epistasisTest(const Snp & other)const
{
    float retVal=0.0;

    ContingencyTable t;
        t.M_.fill(0);
    
    //PLUS ONE FOR CORRECTION
    for(int i=0; i<GENOTYPE_LEVELS; ++i)
    {
        //Count controls
        for(int j=0; j<GENOTYPE_LEVELS; ++j)
        {
            t.M_[i+(GENOTYPE_LEVELS*j)] = popCountAnd(allSamples_, other.allSamples_, i*CONR_, j*CONR_, CONR_)+1;
            
        }
        //Count Cases
        for(int j=0; j<GENOTYPE_LEVELS; ++j)
        {
            t.M_[i+(GENOTYPE_LEVELS*j)+(GENOTYPE_LEVELS*GENOTYPE_LEVELS)] = popCountAnd(allSamples_, other.allSamples_, CASS_+(i*CASR_), CASS_+(j*CASR_), CASR_ )+1;
        }
    }   
    
    t.calculateTotals();

    //18 and 27 offset for row and column totals respectively
    for(int y=0; y< GENOTYPE_PAIRINGS; ++y)
    {

        for(int x=0; x< CONTINGENCY_COLUMNS; ++x)
        {            
            //Use margin values
            int c1 = t.M_[y+18];
            int c2 = t.M_[x+27];
            float c3 = (t.M_[27] + t.M_[28]);


            float c4 = t.M_[y+(9*x)];

            float expected  = ( c1 * c2 ) / c3 ;

            retVal = retVal +  (pow(c4-expected,2)/expected);
           
        }
    }

    return retVal;
}

std::array<std::vector<PACK_TYPE>, GENOTYPE_LEVELS> Snp::packGenotypesHelper(const std::vector<char> & genotypes)
{
    std::array<std::vector<PACK_TYPE>, GENOTYPE_LEVELS> packedGenotypes;

    std::array<PACK_TYPE,GENOTYPE_LEVELS> sectionCode{0,0,0};
    for(int i=0; i<genotypes.size(); ++i)
    {
        
        
  
        PACK_TYPE mask = (PACK_TYPE)1<<(PACK_SIZE-(i+1)%PACK_SIZE);
        
        //Also convert from ascii value to numerical
        switch(genotypes[i])
        {
            case '0':
                sectionCode[0] = sectionCode[0] | mask;
                break;
            case '1':
                sectionCode[1] = sectionCode[1] | mask;
                break;
            case '2':
                sectionCode[2] = sectionCode[2] | mask;
                break;
        }
        
        //When a section is full, push them back and reset the sections to empty;
        if( ((i+1) % PACK_SIZE == 0) || ( (i+1) == genotypes.size()) )
        {
            packedGenotypes[0].push_back(sectionCode[0]);
            packedGenotypes[1].push_back(sectionCode[1]);
            packedGenotypes[2].push_back(sectionCode[2]);
            
            sectionCode = {0,0,0};
        }
    }

    return packedGenotypes;
}

void Snp::packGenotypes(std::vector<PACK_TYPE> & allPackedGenotypes, const std::vector<char> & controlGenotypes, const std::vector<char> & caseGenotypes)
{
    std::array<std::vector<PACK_TYPE>, GENOTYPE_LEVELS> packedControls = packGenotypesHelper(controlGenotypes);
    std::array<std::vector<PACK_TYPE>, GENOTYPE_LEVELS> packedCases = packGenotypesHelper(caseGenotypes);

    for(int coding=0; coding<GENOTYPE_LEVELS; ++coding)
        for(int i=0; i< packedControls[0].size(); ++i)
        {
            allPackedGenotypes.push_back(packedControls[coding][i]);
        }

    for(int coding=0; coding<GENOTYPE_LEVELS; ++coding)
        for(int i=0; i< packedCases[0].size(); ++i)
        {
            allPackedGenotypes.push_back(packedCases[coding][i]);
        }
}

