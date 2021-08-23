#include "Snp.h"

//Define static members, set in main
SnpDimensions Snp::dim;
//ID_Sample Snp::COW_;
//ID_Sample Snp::CAW_;
//ID_Sample Snp::CAS_;

/*Snp::Snp(){
    index_ = -1;
    numControls_ = 0;
    numCases_ = 0;
    allSamples_ = std::vector<PackedGenotype>();
}*/

Snp::Snp(ID_Snp index, const GenotypeMatrix& controls, const GenotypeMatrix& cases){
    packGenotypes(controls.rowBegin(index), controls.rowEnd(index), allSamples_);
    packGenotypes(cases.rowBegin(index), cases.rowEnd(index), allSamples_);

    index_ = index;

    //Determine number of vectorized elements
    //COW_ = CONR_ / 4 + (CONR_ % 4 != 0);
    //CAW_ = CASR_ / 4 + (CASR_ % 4 != 0);
    //CAS_ = (3 * COW_);
}

/*Snp::Snp(const Snp & cpy){
    index_ = cpy.index_;
    allSamples_ = cpy.allSamples_;
}*/

Snp::Snp(const Snp & s1, const Snp & s2){
    //Index is -1 because any node created this way will be internal
    //and not refer to a specific snp from the dataset
    index_ = -1; //TODO: FIX TO BE EXPLICIT TYPE   (currently should wrap arount to largest uval)
  //  std::cout << index_ << " ? "<< std::endl;
    allSamples_ = bitMerge(s1.allSamples_, s2.allSamples_);
}

ID_Snp Snp::getIndex()const{
    return index_;
}

float Snp::computeMinorAlleleFrequency()const{
    ID_Sample homoMajor = popCount(allSamples_, 0, dim.CONR_) + popCount(allSamples_, dim.CASS_, dim.CASR_) ;
    ID_Sample hetero = popCount(allSamples_, dim.CONR_, dim.CONR_) + popCount(allSamples_, dim.CASS_+ dim.CASR_, dim.CASR_);
    ID_Sample homoMinor = popCount(allSamples_, 2* dim.CONR_, dim.CONR_) + popCount(allSamples_, dim.CASS_+(2* dim.CASR_), dim.CASR_) ;
    
    return (2*homoMinor+hetero)/(float)(2*homoMajor + hetero + 2*homoMinor);
}

ID_Sample Snp::computeDifferences(const Snp & other)const{
    return ((dim.numControls_ + dim.numCases_)- popCountAnd(allSamples_, other.allSamples_, 0, allSamples_.size()));
}

float Snp::computeUnknownRatio()const{
    int numberKnown = popCount(allSamples_, 0, allSamples_.size());

    return ( numberKnown/(float)(dim.numControls_ + dim.numCases_) );
}

void Snp::fillTable(CTable2& t, const Snp& snp1, const Snp& snp2){
    t.zero();

    for (const auto& e : CTable::rowOrder) {
        const uint8_t rowI = 9 * e.first + e.second * 3;

        t.data[rowI] = popCountAnd_it(snp1.allSamples_.begin() + (e.first * dim.CONR_), snp2.allSamples_.begin() + (e.second * dim.CONR_), dim.CONR_) + 1;
        t.data[rowI + 1] = popCountAnd_it(snp1.allSamples_.begin() + (e.first * dim.CASR_ + dim.CASS_), snp2.allSamples_.begin() + (e.second * dim.CASR_ + dim.CASS_), dim.CASR_) + 1;

        //t.data[rowI] =   popCountAnd_Vec(samples_, other.samples_, e.first * COW_, e.second * COW_, COW_ )+ 1;
       // t.data[rowI + 1] = popCountAnd_Vec(samples_, other.samples_, e.first * CAW_ + CAS_, e.second * CAW_ + CAS_, CAW_) + 1;

        //fill the row totals
        t.data[rowI + 2] = t.data[rowI] + t.data[rowI + 1];

        //accumulate the column totals
        t.data[27] += t.data[rowI];
        t.data[28] += t.data[rowI + 1];
    }
}

float Snp::marginalTest()const{
   float retVal=0.0;
   CTable1 t;
   t.zero();

   for (uint8_t i = 0; i <= 6; i+=3) {
       //Fill cells for this row, correcting by adding 1
       t.data[i] = popCount(allSamples_, i/3 * dim.CONR_, dim.CONR_) + 1;
       t.data[i + 1] = popCount(allSamples_, dim.CASS_ + (i/3 * dim.CASR_), dim.CASR_) + 1;

       //fill the row totals
       t.data[i + 2] = t.data[i] + t.data[i + 1];

       //accumulate the column totals
       t.data[9] += t.data[i];
       t.data[10] += t.data[i + 1];
   }

    constexpr uint8_t factors = 3;
    constexpr uint8_t levels = 2;

    float chi2 = 0.0;
    for (uint8_t row = 0; row < factors; ++row) {
        for (uint8_t col = 0; col < levels; ++col) {
            ID_Sample c1 = t.data[row * 3 +2];
            ID_Sample c2 = t.data[9 + col];
            ID_Sample c3 = t.data[9] + t.data[10];
            ID_Sample c4 = t.data[row * 3 + col];

            float expected  = (c1 * c2) / (float)c3 ;
            chi2 += (pow(c4 - expected, 2) / expected);
        }
    }
    return chi2;
}

void Snp::packGenotypes(std::vector<uint8_t>::const_iterator begin, std::vector<uint8_t>::const_iterator end, std::vector<uint64_t>& dest ){
    std::array<std::vector<PackedGenotype>,3> packedGenotypes;
    std::array<PackedGenotype, 3> sectionCode{ 0,0,0 };

    PackedGenotype offset=0;

    for (auto it = begin; it != end; ++it) {
        PackedGenotype mask = (PackedGenotype)1 << (PACKED_SIZE - (offset + 1) % PACKED_SIZE);

        sectionCode[*it] = sectionCode[*it] | mask;
        
        //When a section is full, push them back and reset the sections to empty;
        if (((offset + 1) % PACKED_SIZE == 0)) {
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


    //pad for AVX
    /*
    while (packedGenotypes[0].size() % 4 != 0) {
        packedGenotypes[0].push_back(0);
        packedGenotypes[1].push_back(0);
        packedGenotypes[2].push_back(0);
    }
    */
    /*
    for (size_t i = 0; i < packedGenotypes[0].size(); i+=4) {
        samples_.push_back(_mm256_loadu_si256((__m256i*) & packedGenotypes[0][i]));
    }
    for (size_t i = 0; i < packedGenotypes[0].size(); i += 4) {
        samples_.push_back(_mm256_loadu_si256((__m256i*) & packedGenotypes[1][i]));
    }
    for (size_t i = 0; i < packedGenotypes[0].size(); i += 4) {
        samples_.push_back(_mm256_loadu_si256((__m256i*) & packedGenotypes[2][i]));
    }
    */
    //------





  
    
}

void Snp::setDimensions(ID_Sample controls, ID_Sample cases) {
    Snp::dim = {
       controls,
       cases,
       //Determine the number of packed elements required to store num samples, adds one element to handle the last non-full element
       controls / PACKED_SIZE + (controls % PACKED_SIZE != 0),
       cases / PACKED_SIZE + (cases % PACKED_SIZE != 0),
       3 * (controls / PACKED_SIZE + (controls % PACKED_SIZE != 0))
    };
}


void Snp::to_serial(std::ostream& os, const Snp& snp) {
    os.write(reinterpret_cast<const char*>(&snp.index_),sizeof(ID_Snp));
    os.write(reinterpret_cast<const char*>(&snp.allSamples_), snp.allSamples_.size()*sizeof(PackedGenotype));
}