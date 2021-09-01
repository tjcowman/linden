#include "TopSnpList.h"

TopSnpList::TopSnpList(ID_Snp numberSnps) :
    topK_(numberSnps),
    cutoff_(0.0),
    currentPartners_(std::vector<std::pair<ID_Snp, float>>(numberSnps, { -1,0.0 })),
    testCounter_({0,0}),
    insertSincePrefix_(0){

    pairwiseSignificanceCounts_.fill(0);
}

TopSnpList::TopSnpList(ID_Snp topK, ID_Snp numberSnps, float cutoff)
{
    topK_ = topK;
    cutoff_ = cutoff;

    currentPartners_ = std::vector<std::pair<ID_Snp, float>>(numberSnps,{-1,0.0});

    testCounter_ = { 0,0 };
    
    pairwiseSignificanceCounts_.fill(0);
    insertSincePrefix_ = 0;
}

bool TopSnpList::attemptInsert(ID_Snp snpIndex1, ID_Snp snpIndex2, float score)
{ 
    bool retVal = false;
    #pragma omp critical
    {
        if(insertSincePrefix_ == topK_){
            ID_Snp reverseSum = 0;
            //Performs a prefix sum from higher to lower because want number of pairs higher than index, not the number index higher than
            for(int i=MAX_CUTOFF; i>getCutoff(); --i){
                reverseSum += pairwiseSignificanceCounts_[i];
                if (reverseSum > topK_) {
                    cutoff_ = i;
                    break;
                }
            }
            insertSincePrefix_ = 0;
        }

        //If the new score is already < the cutoff we dont care about counting it or updating currentPartners_
        if (score > getCutoff()) {
            ++insertSincePrefix_;

            int index = std::min((int)score, MAX_CUTOFF);
            ++pairwiseSignificanceCounts_[index];

            if (pairwiseSignificanceCounts_[index] >= topK_ )//|| pairwiseSignificanceCountsPrefixSum_[index] >= topK_)
                cutoff_ = index;

            //Check if the new pair score is greater than the current best pairs for each constituent SNP
            if (score > currentPartners_[snpIndex1].second) {
                currentPartners_[snpIndex1] = { snpIndex2, score };
                retVal = true;
            }

            if (score > currentPartners_[snpIndex2].second) {
                currentPartners_[snpIndex2] = { snpIndex1, score };
                retVal = true;
            }
        }
    }

    return retVal;
}

ID_Snp TopSnpList::size()const {
    return currentPartners_.size();
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
    //Iterate over the currently found top partners and push the cutoff pasing ones
    for(size_t i=0; i< currentPartners_.size(); ++i)
        if(currentPartners_[i].second >= getCutoff()-1 && currentPartners_[i].second >0)
            formattedPairs_.cutoff.push_back(  TopPairing(i, currentPartners_[i].first, currentPartners_[i].second )       );
    
    sort(formattedPairs_.cutoff.begin(), formattedPairs_.cutoff.end(), TopPairing::orderByScore);

    calculateReciprocalPairs();  
    sort(formattedPairs_.reciprocal.begin(), formattedPairs_.reciprocal.end(), TopPairing::orderByScore);
}

const FormattedPairs& TopSnpList::getPairs()const{
    return formattedPairs_;
}

void TopSnpList::calculateReciprocalPairs(){
    for (size_t i = 1; i < formattedPairs_.cutoff.size(); ++i) {
        if (formattedPairs_.cutoff[i] == formattedPairs_.cutoff[i - 1])
            formattedPairs_.reciprocal.push_back(TopPairing::normalize(formattedPairs_.cutoff[i]));
    }
}
