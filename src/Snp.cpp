#include "Snp.h"

//Define static members, will be set correctly in the Snp constructors
ID_Sample Snp::numControls_;
ID_Sample Snp::numCases_;
int Snp::CONR_; //Number of packed control samples
int Snp::CASR_; //Number of packed case samples
int Snp::CASS_; //Offest for the beginning of the case samples

Snp::Snp(ID_Snp index, const GenotypeMatrix& controls, const GenotypeMatrix& cases){
    packGenotypes(controls.rowBegin(index), controls.rowEnd(index), allSamples_);
    packGenotypes(cases.rowBegin(index), cases.rowEnd(index), allSamples_);

    index_ = index;
    numControls_ = controls.width;
    numCases_ = cases.width;
    
    //Determine the number of packed elements required to store num samples, adds one element to handle the last non-full element
    CONR_ = numControls_ / PACKED_SIZE + (numControls_ % PACKED_SIZE != 0); 
    CASR_ = numCases_ / PACKED_SIZE + (numCases_ % PACKED_SIZE != 0);
    
    CASS_ = (3 * CONR_);//+ 1; 
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
    index_ = -1; //TODO: FIX TO BE EXPLICIT TYPE   (currently should wrap arount to largest uval)
    std::cout << index_ << " ? "<< std::endl;
    allSamples_ = bitMerge(s1.allSamples_, s2.allSamples_);
}

ID_Snp Snp::getIndex()const
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


ID_Sample Snp::computeDifferences(const Snp & other)const
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
   CTable1 t;
    
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

    CTable2 t;
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


void Snp::packGenotypes(std::vector<uint8_t>::const_iterator begin, std::vector<uint8_t>::const_iterator end, std::vector<uint64_t>& dest ){
    std::array<std::vector<PackedGenotype>,3> packedGenotypes;
    std::array<PackedGenotype, 3> sectionCode{ 0,0,0 };

    PackedGenotype offset=0;

    for (auto it = begin; it != end; ++it) {
        PackedGenotype mask = (PackedGenotype)1 << (PACKED_SIZE - (offset + 1) % PACKED_SIZE);

        sectionCode[*it] = sectionCode[*it] | mask;
        
        //When a section is full, push them back and reset the sections to empty;
        if (((offset + 1) % PACKED_SIZE == 0))//((offset + 1) == genotypes.size()))
        {
            packedGenotypes[0].push_back(sectionCode[0]);
            packedGenotypes[1].push_back(sectionCode[1]);
            packedGenotypes[2].push_back(sectionCode[2]);

            sectionCode = { 0,0,0 };
            offset = 0;
        }
        else {
            ++offset;
        }
    }

    //Push final sections if not full
    if (sectionCode != std::array<PackedGenotype, 3>{0, 0, 0}) {
        packedGenotypes[0].push_back(sectionCode[0]);
        packedGenotypes[1].push_back(sectionCode[1]);
        packedGenotypes[2].push_back(sectionCode[2]);
    }

    dest.insert(dest.end(), packedGenotypes[0].begin(), packedGenotypes[0].end());
    dest.insert(dest.end(), packedGenotypes[1].begin(), packedGenotypes[1].end());
    dest.insert(dest.end(), packedGenotypes[2].begin(), packedGenotypes[2].end());
}