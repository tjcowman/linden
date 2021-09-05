#include "SnpSet.h"

SnpSet::SnpSet() {

}

SnpSet::SnpSet(const std::vector<Locus>& loci, const GenotypeMatrix& controls, const GenotypeMatrix& cases) :
    sizeUnfiltered(loci.size()),
    dim(SnpDimensions(controls.width, cases.width))
{
    for (ID_Snp i = 0; i < loci.size(); ++i) {
        data.push_back(Snp(i, controls, cases));
        //locations.push_back(loci[i].location);
        
    }

    this->loci = loci;
    dim = SnpDimensions(controls.width, cases.width);
    
}

ID_Snp SnpSet::size()const {
    return data.size();
}

ID_Snp SnpSet::getSizeUnfiltered()const {
    return sizeUnfiltered;
}

const SnpDimensions& SnpSet::getDimensions()const {
    return dim;
}

const std::vector<Snp>& SnpSet::getSnps() {
    return data;
}

//const std::vector<Location>& SnpSet::getLocations() {
   // return locations;
//}

void SnpSet::to_serial(std::ostream& os, const SnpSet& snpSet) {
    os.write(reinterpret_cast<const char*>(&snpSet.sizeUnfiltered), sizeof(ID_Snp));
    os.write(reinterpret_cast<const char*>(&snpSet.dim), sizeof(SnpDimensions));
    vector_to_serialc<Snp, ID_Snp>(os, snpSet.data);
   // vector_to_serial<Location, ID_Snp>(os, snpSet.locations);
    vector_to_serialc<Locus, ID_Snp>(os, snpSet.loci);
}

SnpSet SnpSet::from_serial(std::istream& is) {
    SnpSet e;
    is.read(reinterpret_cast<char*>(&e.sizeUnfiltered), sizeof(ID_Snp));
    is.read(reinterpret_cast<char*>(&e.dim), sizeof(SnpDimensions));
    e.data = vector_from_serialc<Snp, ID_Snp>(is);
    //e.locations = vector_from_serial<Location, ID_Snp>(is);
    e.loci = vector_from_serialc<Locus, ID_Snp>(is);

   // std::cerr << "TMP FROM SET SERIAL" << std::endl;
    return e;
}