
#include "SnpSet.h"

#include "Serializers.hpp"

SnpSet::SnpSet(const std::vector<Locus>& loci, const GenotypeMatrix& controls, const GenotypeMatrix& cases) :
    sizeUnfiltered(loci.size()),
    dim(Snp::Dimensions(controls.width, cases.width)),
    loci(loci)
{
    for (ID_Snp i = 0; i < loci.size(); ++i)
    {
        data.push_back(Snp(i, controls, cases));
    }
}

void SnpSet::to_serial(std::ostream& os, const SnpSet& snpSet)
{
    os.write(reinterpret_cast<const char*>(&snpSet.sizeUnfiltered), sizeof(ID_Snp));
    os.write(reinterpret_cast<const char*>(&snpSet.dim), sizeof(Snp::Dimensions));
    vector_to_serialc<Snp, ID_Snp>(os, snpSet.data);
    vector_to_serialc<Locus, ID_Snp>(os, snpSet.loci);
}

SnpSet SnpSet::from_serial(std::istream& is)
{
    SnpSet e;

    auto Bread = is.tellg();
    is.read(reinterpret_cast<char*>(&e.sizeUnfiltered), sizeof(ID_Snp));

    auto Aread = is.tellg();
    is.read(reinterpret_cast<char*>(&e.dim), sizeof(Snp::Dimensions));

    if(is.fail())
    {
        std::cerr<<"Read fail"<<std::endl;
    }

    e.data = vector_from_serialc<Snp, ID_Snp>(is);

    auto Rread = is.tellg();
    e.loci = vector_from_serialc<Locus, ID_Snp>(is);

    auto Lread = is.tellg();

    return e;
}