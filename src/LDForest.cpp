#include "LDForest.h"

LDForest::LDForest( int topKSnps, int numberControlSamples, int numberCaseSamples, int numberSnps)
{
    topSnpList_ = TopSnpList(topKSnps, numberSnps, 0);
}

void LDForest::insert(const Snp & snp, char chromosome, int basePair)
{
    ldtrees_.push_back(LDTree(snp, chromosome, basePair));
}

size_t LDForest::size()const
{
    return ldtrees_.size();
}

void LDForest::mergeTrees(double maxUnknownFraction, DatasetSizeInfo datasetSizeInfo)
{
    std::clog<<"merging LD trees"<<std::endl;

    uint64_t change = UINT64_MAX;//INT_MAX;
    double unknownFraction = 0.0;

    while(change > 5 || unknownFraction < maxUnknownFraction)
    {
        change = mergeTreeIteration(unknownFraction, datasetSizeInfo);

        std::clog<<"\tremaining: "<<size()<<"          \r"<<std::flush;
        
        unknownFraction += .01;
        if(unknownFraction > maxUnknownFraction)
            unknownFraction = maxUnknownFraction;
    }
    std::clog<<"\n";

}

size_t LDForest::mergeTreeIteration(float unknownFraction, DatasetSizeInfo datasetSizeInfo)
{
    std::vector<LDTree> mergedTrees;
    
    int allowedDifferences = unknownFraction * (datasetSizeInfo.controls_ + datasetSizeInfo.cases_);
    
    size_t beforeSize = ldtrees_.size();

    for(size_t i=0; i<size(); ++i)
    {
        if(!ldtrees_[i].empty())
            for(size_t j=i+1; j < std::min(i+10,size()); ++j)
            {
                if((!ldtrees_[j].empty()) && (ldtrees_[i].size() == ldtrees_[j].size()))
                {
                    if(ldtrees_[i].computeDifferences(ldtrees_[j]) < allowedDifferences)
                    {
                        LDTree potentialMerge(ldtrees_[i], ldtrees_[j]);
                        {
                            mergedTrees.push_back(potentialMerge);
                            ldtrees_[i].clear();
                            ldtrees_[j].clear();
                            break;
                        }
                    }
                }
            }
            
        if(!ldtrees_[i].empty())
            mergedTrees.push_back(ldtrees_[i]);
    }
    
    ldtrees_ = mergedTrees;
		
	return beforeSize - ldtrees_.size();
}

void LDForest::testTrees(int maxThreadUsage)
{
    std::clog<<"testing Trees"<<std::endl;
    std::clog<<"\tcompleted: 0/"<<size()<<"               \r"<<std::flush;
    
    size_t treesFinished = 0;
    #pragma omp parallel for num_threads(maxThreadUsage) schedule(dynamic, 20)
    for(int i=0; i<size(); ++i)
    {
		
        //if(parameterInfo.noTrees_ == 0)
        {
            for(int j=i+1; j<size(); ++j)
                ldtrees_[i].epistasisTest(ldtrees_[j], topSnpList_);
        }
        /*else
        {
            for(int j=i+1; j<size(); ++j)
                ldtrees_[i].epistasisTestNoTrees(ldtrees_[j], topSnpList_);
        }*/
    
        if(i % 100 == 0)
        {   
            #pragma omp critical
            {
            treesFinished +=100;
            std::clog<<"\tcompleted: "<<std::min(treesFinished, size())<<"/"<<size()<<"               \r"<<std::flush;
            }
                
        }
            
    }
    std::clog<<std::endl;
}

void LDForest::writeResults(const std::vector<Locus>& infoMatrix, Args& args, DatasetSizeInfo datasetSizeInfo)
{
    topSnpList_.calculateFormattedResults();
    
    std::ofstream ofs;
    
    if(args.output != "")
    {
        ofs.open(args.output + ".reciprocalPairs", std::ofstream::out); // | ofstream::app);
        {
            auto rp = topSnpList_.getReciprocalPairs();
            std::sort(rp.begin(), rp.end(), TopPairing::orderByScore);
            ofs << "chi2\tid1\tchr1\tbp1\tid2\tchr2\tbp2\n";
            for(const auto& e : rp)
            {
                ofs << e.score_ << "\t" << infoMatrix[e.indexes_.first] <<"\t"<< infoMatrix[e.indexes_.second]<<"\n";
            }
            
        }
        ofs.close();
        
        ofs.open(args.output + ".cutoffPairs", std::ofstream::out);
        {
            auto rp = topSnpList_.getCutoffPairs();
            std::sort(rp.begin(), rp.end(), TopPairing::orderByScore);
            ofs << "chi2\tid1\tchr1\tbp1\tid2\tchr2\tbp2\n";
            for(const auto& e : rp)
            {
                ofs << e.score_ << "\t" << infoMatrix[e.indexes_.first] << "\t"<< infoMatrix[e.indexes_.second] << "\n";
            }
        }
        ofs.close();
    }
    else
    {
        auto rp = topSnpList_.getReciprocalPairs();
        std::sort(rp.begin(), rp.end(), TopPairing::orderByScore);
        for(const auto& e : rp)
        {
            std::cout << e.score_ << "\t" << infoMatrix[e.indexes_.first] << "\t" << infoMatrix[e.indexes_.second] << "\trecip\n";
        }
        
        rp = topSnpList_.getCutoffPairs();
        std::sort(rp.begin(), rp.end(), TopPairing::orderByScore);
        for(const auto& e : rp)
        {
            std::cout << e.score_ << "\t" << infoMatrix[e.indexes_.first] << "\t" << infoMatrix[e.indexes_.second] << "\tcutoff\n";
        }
        
    }
    
    std::clog<<"finished"<<"\n";
    std::clog<<"\tleaf tests: "<<topSnpList_.getLeafTests()<<"\n";
    std::clog<<"\tinternal tests: "<<topSnpList_.getInternalTests()<<"\n";
}
