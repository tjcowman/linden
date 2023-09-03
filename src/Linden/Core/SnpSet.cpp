
#include "SnpSet.h"

namespace Linden::Core
{
    SnpSet::SnpSet(const std::vector<Genetics::Locus>& loci, const Genetics::GenotypeMatrix& controls, const Genetics::GenotypeMatrix& cases) :
        sizeUnfiltered(loci.size()),
        dim(Snp::Dimensions(controls.width, cases.width)),
        loci(loci)
    {
        for (Genetics::Id::Snp i = 0; i < loci.size(); ++i)
        {
            data.push_back(Snp(i, controls, cases));
        }
    }

    void SnpSet::to_serial(std::ostream& os, const SnpSet& snpSet)
    {
        // TODO Implement
        // os.write(reinterpret_cast<const char*>(&snpSet.sizeUnfiltered), sizeof(Genetics::Id::Snp));
        // os.write(reinterpret_cast<const char*>(&snpSet.dim), sizeof(Snp::Dimensions));
        // vector_to_serialc<Snp, Genetics::Id::Snp>(os, snpSet.data);
        // vector_to_serialc<Genetics::Locus, Genetics::Id::Snp>(os, snpSet.loci);
    }

    SnpSet SnpSet::from_serial(std::istream& is)
    {
        // TODO Implement
        SnpSet e;

        /*auto Bread = is.tellg();
        is.read(reinterpret_cast<char*>(&e.sizeUnfiltered), sizeof(Genetics::Id::Snp));

        auto Aread = is.tellg();
        is.read(reinterpret_cast<char*>(&e.dim), sizeof(Snp::Dimensions));

        if(is.fail())
        {
            std::cerr<<"Read fail"<<std::endl;
        }

        e.data = vector_from_serialc<Snp, Genetics::Id::Snp>(is);

        auto Rread = is.tellg();
        e.loci = vector_from_serialc<Genetics::Locus, Genetics::Id::Snp>(is);

        auto Lread = is.tellg();
*/
        return e;
    }
}