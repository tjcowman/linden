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

static const int MAX_CUTOFF = 100;//The maximum chi2 value to differentiate between
//static const int PREFIX_SUM_ROLLOVER = 1000;

#include "Types.h"


//Used in the analysis and output of top pairs, not in discovery and cutoff updating
struct TopPairing{
    TopPairing(ID_Snp snpIndex1, ID_Snp snpIndex2, float score) : indexes_({snpIndex1,snpIndex2}), score_(score) {
        score_ = score;
        indexes_ = { snpIndex1, snpIndex2 };
    }

    /**
    * Creates a new TopPairing such that the smaller index is in indexes_.first for reciprocal detection and consistency
    */
    static TopPairing normalize(const TopPairing& p) {
        if(p.indexes_.first < p.indexes_.second)
            return TopPairing{
                p.indexes_.first,
                p.indexes_.second,
                p.score_,
            };
        else
            return TopPairing{
                p.indexes_.second,
                p.indexes_.first,
                p.score_,
            };
    }
    

    bool operator==(const TopPairing& rhs)const {
        return((indexes_.first == rhs.indexes_.first || indexes_.first == rhs.indexes_.second) &&
            (indexes_.second == rhs.indexes_.first || indexes_.second == rhs.indexes_.second)
            );
    }

    void normalize() {
        if (indexes_.first < indexes_.second)
            std::swap(indexes_.first, indexes_.second);
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

struct FormattedPairs {
    std::vector<TopPairing> cutoff;
    std::vector<TopPairing> reciprocal;
};

class TopSnpList{
    public:
        /**
        * IMPORTANT : Currently needs numberSnps to be set to the size of the possible SNPIndexes(BEFORE ANY FILTERING)
        */
        TopSnpList(ID_Snp topK, ID_Snp numberSnps, float cutoff);          
        bool attemptInsert(ID_Snp snpIndex1, ID_Snp snpIndex2, float score);
        
        float getCutoff()const;
        
        void incrementTestCounter(const TestCounter& tests);    
        const TestCounter& getTestCounter()const;

        void calculateFormattedResults();
        const FormattedPairs& getPairs()const;
 
    private:
        void calculateReciprocalPairs();

        //Stores the most up to date best pairing for each SNP
        std::vector<std::pair<ID_Snp, float>> currentPartners_;

        TestCounter testCounter_;
        
        float cutoff_;
        int topK_;

        /**
        * Holds the count of detected pairs with <index> significance
        */
        std::array<int, MAX_CUTOFF+1> pairwiseSignificanceCounts_;
    
        ID_Snp  insertSincePrefix_;

        //Formatted results
        FormattedPairs formattedPairs_;
    
};
#endif //TOP_SNP_LIST_H
