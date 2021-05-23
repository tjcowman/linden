#include "TopSnpList.h"

TopSnpList::TopSnpList()
{
    
}

TopSnpList::TopSnpList(const TopSnpList & cpy)
{
    topK_ = cpy.topK_;
    cutoff_ = cpy.cutoff_;
    
    topScores_ = cpy.topScores_;
    topPartners_ = cpy.topPartners_;
    internalTestsCounter_ = cpy.internalTestsCounter_;
    leafTestsCounter_ = cpy.leafTestsCounter_;
    
    pairwiseSignificanceCounts_ = cpy.pairwiseSignificanceCounts_;
    pairwiseSignificanceCountsPrefixSum_ = cpy.pairwiseSignificanceCountsPrefixSum_;
    prefixSumTimer_ = cpy.prefixSumTimer_;
}

TopSnpList::TopSnpList(int topK, int numberSnps, float cutoff)
{
    topK_ = topK;
    cutoff_ = cutoff;
    
    topScores_ = std::vector<float>(numberSnps, 0.0);
    topPartners_ = std::vector<int>(numberSnps, -1);
    internalTestsCounter_ = 0;
    leafTestsCounter_ = 0;
    
    pairwiseSignificanceCounts_.fill(0);
    pairwiseSignificanceCountsPrefixSum_.fill(0);
    prefixSumTimer_ = 0;
}

bool TopSnpList::attemptInsert(int snpIndex1, int snpIndex2, float score)
{ 
    bool retVal = false;

    #pragma omp critical
    {
        if( prefixSumTimer_ >= PREFIX_SUM_ROLLOVER )
        {
            pairwiseSignificanceCountsPrefixSum_ = pairwiseSignificanceCounts_;
            for(int i=MAX_CUTOFF; i>getCutoff(); --i)
            {
                pairwiseSignificanceCountsPrefixSum_[i-1] = pairwiseSignificanceCountsPrefixSum_[i-1]+pairwiseSignificanceCounts_[i];
            }
            
            prefixSumTimer_ = 0;
        }
        
        if(score > getCutoff())
        {
            int index = std::min((int)score, MAX_CUTOFF);
            pairwiseSignificanceCounts_[index]++;
            
            if(pairwiseSignificanceCounts_[index] >= topK_ || pairwiseSignificanceCountsPrefixSum_[index] >= topK_)
                cutoff_ = index;
        } 
        
        if(score > topScores_[snpIndex1])
        {
            topScores_[snpIndex1] = score;
            topPartners_[snpIndex1] = snpIndex2;
            retVal = true;
        }
        
        if(score > topScores_[snpIndex2])
        {
            topScores_[snpIndex2] = score;
            topPartners_[snpIndex2] = snpIndex1;
            retVal = true;
        }
        prefixSumTimer_++;
    }

    return retVal;
}

float TopSnpList::getCutoff()const
{
    return cutoff_;
}

void TopSnpList::incrementInternalTestsCounter(uint64_t testsDone)
{
    #pragma omp critical
    internalTestsCounter_ += testsDone;
}

void TopSnpList::incrementLeafTestsCounter(uint64_t testsDone)
{
    #pragma omp critical
    leafTestsCounter_ += testsDone;
}

uint64_t TopSnpList::getInternalTests()const
{
    return internalTestsCounter_;
}

uint64_t TopSnpList::getLeafTests()const
{
    return leafTestsCounter_;
}

void TopSnpList::calculateFormattedResults()
{
    
    for(int i=0; i< topPartners_.size(); ++i)
        if(topScores_[i] >= getCutoff()-1 && topScores_[i] >0)
            cutoffPairs_.push_back(  TopPairing(i, topPartners_[i], topScores_[i] )       );

    
    calculateReciprocalPairs();
    
    sort(reciprocalPairs_.begin(), reciprocalPairs_.end(), TopPairing::orderByScore);
    sort(cutoffPairs_.begin(), cutoffPairs_.end(), TopPairing::orderByScore);
}

const std::vector<TopPairing>& TopSnpList::getReciprocalPairs()const
{
    return reciprocalPairs_;
}

const std::vector<TopPairing>& TopSnpList::getCutoffPairs()const
{
    return cutoffPairs_;
}

int TopSnpList::getNumberOfReciprocalPairs()const
{
    return reciprocalPairs_.size();
}

std::ostream& operator<< (std::ostream &out, const TopSnpList & topSnpList)
{
    std::vector< std::pair<float, std::pair<int,int> > > passingPairs;
    
    for(int i=0; i< topSnpList.topPartners_.size(); ++i)
    {
        
        if(topSnpList.topScores_[i] >= topSnpList.getCutoff()-1)
        {
            //cout<< topSnpList.topScores_[i]<<endl;
            passingPairs.push_back(std::make_pair( topSnpList.topScores_[i], std::make_pair(i, topSnpList.topPartners_[i])       )       );
        }
        
    }
    
    sort(passingPairs.rbegin(), passingPairs.rend());

    std::map<TopPairing, int> reciprocalCheck;
    
    for(int i=0; i< passingPairs.size(); ++i)
    {   
        TopPairing tempPair(passingPairs[i].second.first, passingPairs[i].second.second, passingPairs[i].first);
        
        if(reciprocalCheck.count(tempPair)>0)
            reciprocalCheck[tempPair]++;
        else
        {
            reciprocalCheck.insert(std::make_pair(tempPair, 1));
        } 
    }
    
    std::vector<float> lazySort;
    int recips =0;
    for(auto it=reciprocalCheck.begin(); it!=reciprocalCheck.end(); ++it)
    {
        if(it->second == 2)
        {
            
            recips++;
            //cout<<it->first.score_<<endl;
            lazySort.push_back(it->first.score_);
        }
    }
    
    sort(lazySort.rbegin(), lazySort.rend());
    for(int i=0; i< lazySort.size(); ++i)
        std::cout<<lazySort[i]<<"\n";
    
    
    out<<"INTERNAL TESTS:   "<<topSnpList.internalTestsCounter_<<"\n";
    out<<"LEAF TESTS:       "<<topSnpList.leafTestsCounter_<<"\n";
    out<<"PASSING:          "<<passingPairs.size()<<"\n";
    out<<"RECIPROCALS:      "<<recips<<"\n";
    
    return out;
}

void TopSnpList::calculateReciprocalPairs()
{
    std::map<TopPairing, int> reciprocalCheck;
    for(int i=0; i< topPartners_.size(); ++i)
    {
        if(topScores_[i] >= getCutoff()-1)
        {
            TopPairing tempPair(i, topPartners_[i], topScores_[i]);
            
            if(reciprocalCheck.count(tempPair)>0)
                reciprocalCheck[tempPair]++;
            else
                reciprocalCheck.insert(std::make_pair(tempPair, 1));
        }
    }
    
    for(auto it=reciprocalCheck.begin(); it!=reciprocalCheck.end(); ++it)
    {
        if(it->second == 2)
            reciprocalPairs_.push_back(it->first);
    }
}
