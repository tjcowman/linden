#include "TopSnpList.h"

TopSnpList::TopSnpList(ID_Snp topK, ID_Snp numberSnps, float cutoff)
{
    topK_ = topK;
    cutoff_ = cutoff;

    currentPartners_ = std::vector<std::pair<ID_Snp, float>>(numberSnps,{-1,0.0});

    testCounter_ = { 0,0 };
    
    pairwiseSignificanceCounts_.fill(0);
    pairwiseSignificanceCountsPrefixSum_.fill(0);
    prefixSumTimer_ = 0;
}

bool TopSnpList::attemptInsert(ID_Snp snpIndex1, ID_Snp snpIndex2, float score)
{ 
    bool retVal = false;

    #pragma omp critical
    {
        if( prefixSumTimer_ >= PREFIX_SUM_ROLLOVER ){
            pairwiseSignificanceCountsPrefixSum_ = pairwiseSignificanceCounts_;
            for(int i=MAX_CUTOFF; i>getCutoff(); --i)
            {
                pairwiseSignificanceCountsPrefixSum_[i-1] = pairwiseSignificanceCountsPrefixSum_[i-1]+pairwiseSignificanceCounts_[i];
            }
            
            prefixSumTimer_ = 0;
        }
        
        if(score > getCutoff()) {
            int index = std::min((int)score, MAX_CUTOFF);
            pairwiseSignificanceCounts_[index]++;
            
            if(pairwiseSignificanceCounts_[index] >= topK_ || pairwiseSignificanceCountsPrefixSum_[index] >= topK_)
                cutoff_ = index;
        } 
        
        if(score > currentPartners_[snpIndex1].second){
            currentPartners_[snpIndex1] = {snpIndex2, score};
            retVal = true;
        }
        
        if(score > currentPartners_[snpIndex2].second) {
            currentPartners_[snpIndex2] = { snpIndex1, score };
            retVal = true;
        }
        prefixSumTimer_++;
    }

    return retVal;
}

float TopSnpList::getCutoff()const{
    return cutoff_;
}

void TopSnpList::incrementTestCounter(const TestCounter& tests){
    #pragma omp critical
    {
        testCounter_.internal += tests.internal;
        testCounter_.leaf += tests.leaf;
    }
}

const TestCounter& TopSnpList::getTestCounter()const{
    return testCounter_;
}

void TopSnpList::calculateFormattedResults(){
    for(size_t i=0; i< currentPartners_.size(); ++i)
        if(currentPartners_[i].second >= getCutoff()-1 && currentPartners_[i].second >0)
            cutoffPairs_.push_back(  TopPairing(i, currentPartners_[i].first, currentPartners_[i].second )       );

    
    calculateReciprocalPairs();
    
    sort(reciprocalPairs_.begin(), reciprocalPairs_.end(), TopPairing::orderByScore);
    sort(cutoffPairs_.begin(), cutoffPairs_.end(), TopPairing::orderByScore);
}

const std::vector<TopPairing>& TopSnpList::getReciprocalPairs()const{
    return reciprocalPairs_;
}

const std::vector<TopPairing>& TopSnpList::getCutoffPairs()const{
    return cutoffPairs_;
}

std::ostream& operator<< (std::ostream &out, const TopSnpList & topSnpList){
    std::vector< std::pair<float, std::pair<int,int> > > passingPairs;
    
    for(size_t i=0; i< topSnpList.currentPartners_.size(); ++i) {
        if(topSnpList.currentPartners_[i].second >= topSnpList.getCutoff()-1) {
            passingPairs.push_back(std::make_pair( topSnpList.currentPartners_[i].first, std::make_pair(i, topSnpList.currentPartners_[i].second)       )       );
        }      
    }
    
    sort(passingPairs.rbegin(), passingPairs.rend());

    std::map<TopPairing, int> reciprocalCheck;
    
    for(size_t i=0; i< passingPairs.size(); ++i) {   
        TopPairing tempPair(passingPairs[i].second.first, passingPairs[i].second.second, passingPairs[i].first);
        
        if(reciprocalCheck.count(tempPair)>0)
            reciprocalCheck[tempPair]++;
        else{
            reciprocalCheck.insert(std::make_pair(tempPair, 1));
        } 
    }
    
    std::vector<float> lazySort;
    int recips =0;
    for(auto it=reciprocalCheck.begin(); it!=reciprocalCheck.end(); ++it) {
        if(it->second == 2){
            
            recips++;
            lazySort.push_back(it->first.score_);
        }
    }
    
    sort(lazySort.rbegin(), lazySort.rend());
    for(size_t i=0; i< lazySort.size(); ++i)
        std::cout<<lazySort[i]<<"\n";
    
    
    out<<"INTERNAL TESTS:   "<<topSnpList.testCounter_.internal<<"\n";
    out<<"LEAF TESTS:       "<<topSnpList.testCounter_.leaf<<"\n";
    out<<"PASSING:          "<<passingPairs.size()<<"\n";
    out<<"RECIPROCALS:      "<<recips<<"\n";
    
    return out;
}

void TopSnpList::calculateReciprocalPairs(){
    std::map<TopPairing, int> reciprocalCheck;
    for(size_t i=0; i< currentPartners_.size(); ++i){
        if(currentPartners_[i].second >= getCutoff()-1){
            TopPairing tempPair(i, currentPartners_[i].first, currentPartners_[i].second);
            
            if(reciprocalCheck.count(tempPair)>0)
                reciprocalCheck[tempPair]++;
            else
                reciprocalCheck.insert(std::make_pair(tempPair, 1));
        }
    }
    
    for(auto it=reciprocalCheck.begin(); it!=reciprocalCheck.end(); ++it){
        if(it->second == 2)
            reciprocalPairs_.push_back(it->first);
    }
}
