#include "Snp.hpp"

#include "Serializers.hpp"

//Define static members, set in main
SnpDimensions Snp::dim;

// Construct a Snp from an index and vector of samples, ex: from serialized data
Snp::Snp(ID_Snp index, const std::vector<Bitwise::Genotype>& samples) :
    index_(index),
    allSamples_(samples)
{ }

////////////////////////////////////////////////////////////////////////
//! Constructs a new Snp with the specifed index and a matrix for
//! control and case genotypes. Note the index is used to look up the
//! genotype rows in the matrices.
//!
//! @param index The row in the genotype matrices.
//! @param controls Control sample genotype matrix.
//! @param cases Case sample genotype matrix.
////////////////////////////////////////////////////////////////////////
Snp::Snp(ID_Snp index, const GenotypeMatrix& controls, const GenotypeMatrix& cases)
{
    packGenotypes(controls.rowBegin(index), controls.rowEnd(index), allSamples_);
    packGenotypes(cases.rowBegin(index), cases.rowEnd(index), allSamples_);
    index_ = index;
}

//TODO: FIX TO BE EXPLICIT TYPE (currently should wrap arount to largest uval)
////////////////////////////////////////////////////////////////////////
//! Constructs a new snp by merging two children snps. Note that the
//! index in this case is -1 as it will be an internal snp.
//!
//! @param s1 The left snp node.
//! @param s2 The right snp node.
////////////////////////////////////////////////////////////////////////
Snp::Snp(const Snp & s1, const Snp & s2) :
    index_(-1),
    allSamples_(Bitwise::merge(s1.allSamples_, s2.allSamples_))
{ }

float Snp::computeMinorAlleleFrequency() const
{
    ID_Sample homoMajor = Bitwise::count(allSamples_.begin(), dim.CONR_) + Bitwise::count(allSamples_.begin() + dim.CASS_, dim.CASR_) ;
    ID_Sample hetero = Bitwise::count(allSamples_.begin() + dim.CONR_, dim.CONR_) + Bitwise::count(allSamples_.begin() + dim.CASS_+ dim.CASR_, dim.CASR_);
    ID_Sample homoMinor = Bitwise::count(allSamples_.begin() + (2* dim.CONR_), dim.CONR_) + Bitwise::count(allSamples_.begin() + dim.CASS_+(2* dim.CASR_), dim.CASR_) ;

    return (2*homoMinor+hetero)/(float)(2*homoMajor + hetero + 2*homoMinor);
}

ID_Sample Snp::computeDifferences(const Snp & other) const
{
    return ((dim.numControls_ + dim.numCases_)- Bitwise::andCount(allSamples_.begin(), other.allSamples_.begin(), allSamples_.size()));
}

float Snp::computeUnknownRatio() const
{
    int numberKnown = Bitwise::count(allSamples_.begin(), allSamples_.size());

    return ( numberKnown/(float)(dim.numControls_ + dim.numCases_) );
}

void Snp::fillTable(ContingencyTable2& t, const Snp& snp1, const Snp& snp2)
{
    t.zero();

    for (const auto& e : ContingencyTable::rowOrder) {
        const uint8_t rowI = 9 * e.first + e.second * 3;

        t.data[rowI] = Bitwise::andCount(snp1.allSamples_.begin() + (e.first * dim.CONR_), snp2.allSamples_.begin() + (e.second * dim.CONR_), dim.CONR_) + 1;
        t.data[rowI + 1] = Bitwise::andCount(snp1.allSamples_.begin() + (e.first * dim.CASR_ + dim.CASS_), snp2.allSamples_.begin() + (e.second * dim.CASR_ + dim.CASS_), dim.CASR_) + 1;

        //fill the row totals
        t.data[rowI + 2] = t.data[rowI] + t.data[rowI + 1];

        //accumulate the column totals
        t.data[27] += t.data[rowI];
        t.data[28] += t.data[rowI + 1];
    }
}

float Snp::marginalTest() const
{
   ContingencyTable1 t;
   t.zero();

   for (uint8_t i = 0; i <= 6; i+=3) {
       //Fill cells for this row, correcting by adding 1
       t.data[i] = Bitwise::count(allSamples_.begin() + (i/3 * dim.CONR_), dim.CONR_) + 1;
       t.data[i + 1] = Bitwise::count(allSamples_.begin() + (dim.CASS_ + (i/3 * dim.CASR_)), dim.CASR_) + 1;

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

void Snp::packGenotypes(std::vector<uint8_t>::const_iterator begin, std::vector<uint8_t>::const_iterator end, std::vector<Bitwise::Genotype>& dest )
{
    std::array<std::vector<Bitwise::Genotype>,3> packedGenotypes;
    std::array<Bitwise::Genotype, 3> sectionCode{ 0,0,0 };

    Bitwise::Genotype offset=0;

    for (auto it = begin; it != end; ++it) {
        Bitwise::Genotype mask = (Bitwise::Genotype)1 << (Bitwise::Size - (offset + 1) % Bitwise::Size);

        sectionCode[*it] = sectionCode[*it] | mask;

        //When a section is full, push them back and reset the sections to empty;
        if (((offset + 1) % Bitwise::Size == 0)) {
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

    if (sectionCode != std::array<Bitwise::Genotype, 3>{0, 0, 0}) {
        packedGenotypes[0].push_back(sectionCode[0]);
        packedGenotypes[1].push_back(sectionCode[1]);
        packedGenotypes[2].push_back(sectionCode[2]);
    }


    dest.insert(dest.end(), packedGenotypes[0].begin(), packedGenotypes[0].end());
    dest.insert(dest.end(), packedGenotypes[1].begin(), packedGenotypes[1].end());
    dest.insert(dest.end(), packedGenotypes[2].begin(), packedGenotypes[2].end());
}

void Snp::to_serial(std::ostream& os, const Snp& snp)
{
    os.write(reinterpret_cast<const char*>(&snp.index_),sizeof(ID_Snp));
    vector_to_serial<Bitwise::Genotype,ID_Snp>(os, snp.allSamples_);
}

Snp Snp::from_serial(std::istream& is)
{
    auto t1 = value_from_serial<ID_Snp>(is);
    auto t2 =  vector_from_serial<Bitwise::Genotype, ID_Snp>(is);

    return Snp(
        std::move(t1),
        std::move(t2)
    );
}
