
#include "Snp.hpp"

namespace Linden::Core
{

    // Define static members, set in main
    Snp::Dimensions Snp::dim;

    ////////////////////////////////////////////////////////////////////////////////
    //! Construct a Snp from an index and vector of samples.
    //! For example, from serialized data.
    //!
    //! @param index Index for this snp.
    //! @param samples Vector of genotype samples
    ////////////////////////////////////////////////////////////////////////////////
    Snp::Snp(Genetics::Id::Snp index, const std::vector<Bitwise::Genotype>& samples) :
        index_(index),
        allSamples_(samples)
    { }

    ////////////////////////////////////////////////////////////////////////////////
    //! Constructs a new Snp with the specifed index and a matrix for
    //! control and case genotypes. Note the index is used to look up the
    //! genotype rows in the matrices.
    //!
    //! @param index The row in the genotype matrices.
    //! @param controls Control sample genotype matrix.
    //! @param cases Case sample genotype matrix.
    ////////////////////////////////////////////////////////////////////////////////
    Snp::Snp(Genetics::Id::Snp index, const Genetics::GenotypeMatrix& controls, const Genetics::GenotypeMatrix& cases)
    {
        packGenotypes(controls.rowBegin(index), controls.rowEnd(index), allSamples_);
        packGenotypes(cases.rowBegin(index), cases.rowEnd(index), allSamples_);
        index_ = index;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! Constructs a new snp by merging two children snps. Note that the
    //! index in this case is -1 as it will be an internal snp.
    //!
    //! @param s1 The left snp node.
    //! @param s2 The right snp node.
    ////////////////////////////////////////////////////////////////////////////////
    Snp::Snp(const Snp & s1, const Snp & s2) :
        //TODO: FIX TO BE EXPLICIT TYPE (currently should wrap arount to largest uval)
        index_(-1),
        allSamples_(Bitwise::merge(s1.allSamples_, s2.allSamples_))
    { }

    ////////////////////////////////////////////////////////////////////////////////
    //! Computes the minor allele frequency for the Snp.
    //!
    //! @returns The minor allele frequency.
    ////////////////////////////////////////////////////////////////////////////////
    float Snp::computeMinorAlleleFrequency() const
    {
        Genetics::Id::Sample homoMajor = 
            Bitwise::count(GetControlsBegin(0), dim.numPackedControls_) + 
            Bitwise::count(GetCasesBegin(0), dim.numPackedCases_) ;
        Genetics::Id::Sample hetero =
            Bitwise::count(GetControlsBegin(1), dim.numPackedControls_) +
            Bitwise::count(GetCasesBegin(1), dim.numPackedCases_);
        Genetics::Id::Sample homoMinor =
            Bitwise::count(GetControlsBegin(2), dim.numPackedControls_) +
            Bitwise::count(GetCasesBegin(2), dim.numPackedCases_) ;

        return (2*homoMinor+hetero)/(float)(2*homoMajor + hetero + 2*homoMinor);
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! Computes the number of samples different between the snps.
    //!
    //! @param other The snp to compare genotypes to.
    //! @returns The number samples with different genotypes.
    ////////////////////////////////////////////////////////////////////////////////
    Genetics::Id::Sample Snp::computeDifferences(const Snp & other) const
    {
        auto matching = Bitwise::andCount(allSamples_.begin(),
                                        other.allSamples_.begin(),
                                        allSamples_.size());

        return dim.numControls_ + dim.numCases_ - matching;
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! Computes the fraction of non-ambiguous genotypes.
    //!
    //! @returns The fraction of non-ambiguous genotypes.
    ////////////////////////////////////////////////////////////////////////////////
    float Snp::computeUnknownRatio() const
    {
        int numberKnown = Bitwise::count(allSamples_.begin(), allSamples_.size());

        return (numberKnown/(float)(dim.numControls_ + dim.numCases_));
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! Fills out a contingency table with the pairwise genotypes of two Snps.
    //!
    //! @param t The contingency table to fill out.
    //! @param snp1 The first Snp.
    //! @param snp2 The second Snp.
    ////////////////////////////////////////////////////////////////////////////////
    void Snp::fillTable(Statistics::ContingencyTable<9>& t, const Snp& snp1, const Snp& snp2)
    {
        t.zero();
        auto& data = t.Data();

        for (const auto& e : Statistics::ContingencyTableTmp::rowOrder)
        {
            const uint8_t rowI = 9 * e.first + e.second * 3;

            data[rowI] = Bitwise::andCount(
                snp1.GetControlsBegin(e.first),
                snp2.GetControlsBegin(e.second), 
                dim.numPackedControls_) + 1;

            data[rowI + 1] = Bitwise::andCount(
                snp1.GetCasesBegin(e.first),
                snp2.GetCasesBegin(e.second), 
                dim.numPackedCases_) + 1;

            //fill the row totals
            t.RowTotal( 3 * e.first + e.second) = data[rowI] + data[rowI + 1];

            //accumulate the column totals
            t.ColTotal(0) += data[rowI];
            t.ColTotal(1) += data[rowI + 1];
        }
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! Computes the marginal significance of the Snp.
    //!
    //! @returns The marginal signifcance.
    ////////////////////////////////////////////////////////////////////////////////
    float Snp::marginalTest() const
    {
        Statistics::ContingencyTable<3> t;
        auto& data = t.Data();
        t.zero();

        for (uint8_t i = 0; i <= 6; i+=3)
        {
            //Fill cells for this row, correcting by adding 1
            data[i] = Bitwise::count(GetControlsBegin(i/3), dim.numPackedControls_) + 1;
            data[i + 1] = Bitwise::count(GetCasesBegin(i/3), dim.numPackedCases_) + 1;

            //fill the row totals
            t.RowTotal(i/3) = data[i] + data[i + 1];

            //accumulate the column totals
            t.ColTotal(0) += data[i];
            t.ColTotal(1) += data[i + 1];
        }
        return t.Chi2();
    }

    ////////////////////////////////////////////////////////////////////////////////
    //! Packs a range of single byte genotype values into the compressed form.
    //!
    //! @param begin The beginning iterator for the range.
    //! @param end The end iterator for the range.
    //! @param dest The destination iterator for the packed genotype.
    ////////////////////////////////////////////////////////////////////////////////
    void Snp::packGenotypes(std::vector<uint8_t>::const_iterator begin,
                            std::vector<uint8_t>::const_iterator end,
                            std::vector<Bitwise::Genotype>& dest)
    {
        std::array<std::vector<Bitwise::Genotype>,3> packedGenotypes;
        std::array<Bitwise::Genotype, 3> sectionCode{0, 0, 0};

        Bitwise::Genotype offset=0;

        for (auto it = begin; it != end; ++it)
        {
            Bitwise::Genotype mask = (Bitwise::Genotype)1 << (Bitwise::Size - (offset + 1) % Bitwise::Size);

            sectionCode[*it] = sectionCode[*it] | mask;

            //When a section is full, push them back and reset the sections to empty;
            if ((offset + 1) % Bitwise::Size == 0)
            {
                packedGenotypes[0].push_back(sectionCode[0]);
                packedGenotypes[1].push_back(sectionCode[1]);
                packedGenotypes[2].push_back(sectionCode[2]);

                sectionCode = {0, 0, 0};
                offset = 0;
            }
            else
            {
                ++offset;
            }
        }

        //Push final sections if not full
        if (sectionCode != std::array<Bitwise::Genotype, 3>{0, 0, 0})
        {
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
        // TODO Implement
        // os.write(reinterpret_cast<const char*>(&snp.index_),sizeof(Genetics::Id::Snp));
        // vector_to_serial<Bitwise::Genotype,Genetics::Id::Snp>(os, snp.allSamples_);
    }

    Snp Snp::from_serial(std::istream& is)
    {
        // TODO Implement
        // auto t1 = value_from_serial<Genetics::Id::Snp>(is);
        // auto t2 =  vector_from_serial<Bitwise::Genotype, Genetics::Id::Snp>(is);

        // return Snp( std::move(t1), std::move(t2));
        return Snp();
    }
}
