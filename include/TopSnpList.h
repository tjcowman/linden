/**
 * @author Tyler Cowman
 * Stores the current most significant interaction for each SNP. 
 */

#ifndef TOP_SNP_LIST_H
#define TOP_SNP_LIST_H

#include <map>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>

#define MAX_CUTOFF 100
#define PREFIX_SUM_ROLLOVER 1000

#include "Types.h"


//Used in the analysis and output of top pairs, not in discovery and cutoff updating
struct TopPairing{
    TopPairing(ID_Snp snpIndex1, ID_Snp snpIndex2, float score){
        score_ = score;
     
        indexes_ = snpIndex1 < snpIndex2 ?
            std::pair<ID_Snp, ID_Snp>{snpIndex1, snpIndex2}:
            std::pair<ID_Snp, ID_Snp>{snpIndex2, snpIndex1};
    }
      
    bool operator <(const TopPairing &other )const{
        return indexes_<other.indexes_;
    }
    
    static bool orderByScore(const TopPairing & a, const TopPairing &  b){
        return a.score_ > b.score_;
    }
    
    friend std::ostream& operator<< (std::ostream &out, const TopPairing & topPairing){
        out<<topPairing.score_<<" "<<topPairing.indexes_.first<<" "<<topPairing.indexes_.second;
        
        return out;
    }
    
    float score_;
    std::pair<ID_Snp, ID_Snp> indexes_;
};

struct TestCounter {
    uint64_t internal;
    uint64_t leaf;
};

class TopSnpList{
    public:
        //IMPORTANT: Currently needs numberSnps to be set to the size of the possible SNPIndexes (BEFORE ANY FILTERING)
        TopSnpList(ID_Snp topK, ID_Snp numberSnps, float cutoff);
           
        bool attemptInsert(ID_Snp snpIndex1, ID_Snp snpIndex2, float score);
        
        float getCutoff()const;
        

        void incrementTestCounter(const TestCounter& tests);    

        const TestCounter& getTestCounter()const;

        void calculateFormattedResults();

        
        const std::vector<TopPairing>& getReciprocalPairs()const;
        const std::vector<TopPairing>& getCutoffPairs()const;
        
        
        friend std::ostream& operator<< (std::ostream &out, const TopSnpList & topSnpList);
        
    private:
        void calculateReciprocalPairs();

        //Stores the most up to date best pairing for each SNP
        std::vector<std::pair<ID_Snp, float>> currentPartners_;


        TestCounter testCounter_;
        
        float cutoff_;
        int topK_;

        std::array<int, MAX_CUTOFF+1> pairwiseSignificanceCounts_;
        std::array<int, MAX_CUTOFF+1> pairwiseSignificanceCountsPrefixSum_;
        int prefixSumTimer_;
    
        //Formatted results
        std::vector<TopPairing> reciprocalPairs_;
        std::vector<TopPairing> cutoffPairs_;
    
};
#endif //TOP_SNP_LIST_H
