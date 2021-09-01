#include "SnpSet.h"

SnpSet::SnpSet(const std::vector<Locus>& loci, const GenotypeMatrix& controls, const GenotypeMatrix& cases) :
    sizeUnfiltered(loci.size()),
    dim(SnpDimensions(controls.width, cases.width))
{
    for (ID_Snp i = 0; i < loci.size(); ++i) {
        data.push_back(Snp(i, controls, cases));
        locations.push_back(loci[i].location);
    }
    
}


const std::vector<Snp>& SnpSet::getSnps() {
    return data;
}

const std::vector<Location>& SnpSet::getLocations() {
    return locations;
}