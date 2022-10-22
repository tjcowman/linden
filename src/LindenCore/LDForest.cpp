#include "LDForest.h"

LDForest::LDForest(ID_Snp numSnps) :
    topSnpList_(TopSnpList(numSnps)){
}

LDForest::LDForest(SnpSet& snpSet, ID_Snp numSnps) :
    topSnpList_(TopSnpList(numSnps)) {

    const auto& s = snpSet.getSnps();
    const auto& l = snpSet.getLoci();

    ldtrees_.reserve(s.size());

    for (ID_Snp i = 0; i < s.size(); ++i)
    {
        ldtrees_.push_back(LDTree(s[i], l[s[i].getIndex()].location, &topSnpList_));
    }
}

size_t LDForest::size() const{
    return ldtrees_.size();
}

bool LDForest::operator==(const LDForest& lhs) const {
    return ldtrees_ == lhs.ldtrees_;
}

void LDForest::mergeTrees(double maxUnknownFraction){
    std::clog<<"merging LD trees"<<std::endl;

    double unknownFraction = 0.0;
    size_t change;
    do {
     //   std::cout << "TMP MERGE ITER" << std::endl;
        change = mergeTreeIteration(unknownFraction);

        std::clog << "\tremaining: " << size() << "          \r" << std::flush;

        unknownFraction += .01;
        unknownFraction = std::min(unknownFraction, maxUnknownFraction);

    } while (change > 5 || unknownFraction < maxUnknownFraction);

    std::clog<<"\n";
}

size_t LDForest::mergeTreeIteration(float unknownFraction){
    std::vector<LDTree> mergedTrees;
      //mergedTrees.reserve(ldtrees_);

    ID_Snp allowedDifferences = unknownFraction * (Snp::getDimensions().numControls_ + Snp::getDimensions().numCases_); //(log_->controls_ + log_->cases_);
    size_t beforeSize = ldtrees_.size();

    for(size_t i=0; i<size(); ++i){
      //  std::cout << "TMP i " << i << std::endl;
        if (!ldtrees_[i].empty()) { //The ldTree may have been merged from a previous iteration of the i loop
            for (size_t j = i + 1; j < std::min(i + 10, size()); ++j){
                if(ldtrees_[i].validMerge(ldtrees_[j], allowedDifferences)){
               //     std::cout << "TMP TREES " << i <<" " << j << std::endl;
                    LDTree potentialMerge(ldtrees_[i], ldtrees_[j]);
                    mergedTrees.push_back(potentialMerge);
                    ldtrees_[i].clear();
                    ldtrees_[j].clear();
                    break; //If a pair was merged, the i tree is no longer valid for merging
                }
            }
        }

        //If no tree was merged in the j loop with this iteration of i
        if(!ldtrees_[i].empty())
            mergedTrees.push_back(ldtrees_[i]);
    }

   // std::cout << "TMP SWAPPING BUFFER" << std::endl;
    ldtrees_ = mergedTrees;

	return beforeSize - ldtrees_.size();
}

void LDForest::testTrees(int maxThreadUsage){
    std::clog<<"testing Trees"<<std::endl;
    std::clog<<"\tcompleted: 0/"<<size()<<"               \r"<<std::flush;

    size_t treesFinished = 0;
    #pragma omp parallel for num_threads(maxThreadUsage) schedule(dynamic, 20)
    for(std::int32_t i=0; i<size(); ++i){
        for (size_t j = i + 1; j < size(); ++j) {
            ldtrees_[i].epistasisTest(ldtrees_[j]);//, topSnpList_);
        }

        if(i % 100 == 0){
            #pragma omp critical
            {
                treesFinished +=100;
                std::clog<<"\tcompleted: "<<std::min(treesFinished, size())<<"/"<<size()<<"               \r"<<std::flush;
            }
        }
    }
    std::clog<<std::endl;
}

void LDForest::writeResults(const std::vector<Linden::Genetics::Locus>& infoMatrix, const std::string& output){
    topSnpList_.calculateFormattedResults();

    std::ofstream ofs;

    //If no output file was provided
    if(output != ""){
        ofs.open(output + ".reciprocalPairs", std::ofstream::out);
        auto rp = topSnpList_.getPairs().reciprocal;// getReciprocalPairs();
        std::sort(rp.begin(), rp.end(), TopPairing::orderByScore);
        ofs << "chi2\tid1\tchr1\tbp1\tid2\tchr2\tbp2\n";
        for(const auto& e : rp){
            ofs << e.score_ << "\t" << infoMatrix[e.indexes_.first] <<"\t"<< infoMatrix[e.indexes_.second]<<"\n";
        }
        ofs.close();

        ofs.open(output + ".cutoffPairs", std::ofstream::out);
        rp = topSnpList_.getPairs().cutoff;
        std::sort(rp.begin(), rp.end(), TopPairing::orderByScore);
        ofs << "chi2\tid1\tchr1\tbp1\tid2\tchr2\tbp2\n";
        for(const auto& e : rp){
            ofs << e.score_ << "\t" << infoMatrix[e.indexes_.first] << "\t"<< infoMatrix[e.indexes_.second] << "\n";
        }
        ofs.close();
    }
    else{ //print to standard out
        auto rp = topSnpList_.getPairs().reciprocal;
        std::sort(rp.begin(), rp.end(), TopPairing::orderByScore);
        for(const auto& e : rp){
            std::cout << e.score_ << "\t" << infoMatrix[e.indexes_.first] << "\t" << infoMatrix[e.indexes_.second] << "\trecip\n";
        }

        rp = topSnpList_.getPairs().cutoff;
        std::sort(rp.begin(), rp.end(), TopPairing::orderByScore);
        for(const auto& e : rp) {
            std::cout << e.score_ << "\t" << infoMatrix[e.indexes_.first] << "\t" << infoMatrix[e.indexes_.second] << "\tcutoff\n";
        }
    }

    std::clog<<"finished"<<"\n";
    std::clog<<"\tleaf tests: "<<topSnpList_.getTestCounter().leaf<<"\n";
    std::clog<<"\tinternal tests: "<<topSnpList_.getTestCounter().internal<<"\n";
}

void LDForest::to_serial(std::ostream& os, const LDForest& e) {
    //Write the number of Snps to build a new TopSnpList
    ID_Snp numSnps = e.topSnpList_.size();
    os.write(reinterpret_cast<const char*>(&numSnps), sizeof(ID_Snp));

    //IMPORTANT WRITE/READ SET SNP SIZE DIMENSIONS HERE AS THEY ARE STATIC TODO: make this less confusing
    auto dim = Snp::getDimensions();
    os.write(reinterpret_cast<const char*>(&dim), sizeof(Snp::Dimensions));

    //Write number of trees
    ID_Snp numTrees = e.ldtrees_.size();
    os.write(reinterpret_cast<const char*>(&numTrees), sizeof(ID_Snp));
    for(const auto& tree : e.ldtrees_){
        LDTree::to_serial(os, tree);
    }

}

LDForest LDForest::from_serial(std::istream& is)
{
    ID_Snp numSnps;
    Snp::Dimensions dim;
    ID_Snp numTrees;

    is.read(reinterpret_cast<char*>(&numSnps), sizeof(ID_Snp));
    is.read(reinterpret_cast<char*>(&dim), sizeof(Snp::Dimensions));
    is.read(reinterpret_cast<char*>(&numTrees), sizeof(ID_Snp));

    Snp::setDimensions(dim.numControls_, dim.numCases_);
    LDForest e(numSnps);

    e.ldtrees_.reserve(numTrees);
    for (ID_Snp i = 0; i < numTrees; ++i)
        e.ldtrees_.emplace_back(LDTree::from_serial(is));

    return e;

}