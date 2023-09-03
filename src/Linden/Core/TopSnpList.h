/**
 * @author Tyler Cowman
 * Stores the current most significant interaction for each SNP.
 */

#pragma once

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

static const int MAX_CUTOFF = 100;//The maximum chi2 value to differentiate between
//static const int PREFIX_SUM_ROLLOVER = 1000;

#include "Id.hpp"

namespace Linden::Core
{
    //Used in the analysis and output of top pairs, not in discovery and cutoff updating
    struct TopPairing{

        TopPairing(Genetics::Id::Snp snpIndex1, Genetics::Id::Snp snpIndex2, float score) :
            indexes_({snpIndex1,snpIndex2}),
            score_(score)
        { }

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


        bool operator==(const TopPairing& rhs) const {
            return((indexes_.first == rhs.indexes_.first || indexes_.first == rhs.indexes_.second) &&
                (indexes_.second == rhs.indexes_.first || indexes_.second == rhs.indexes_.second)
                );
        }

        void normalize() {
            if (indexes_.first < indexes_.second)
                std::swap(indexes_.first, indexes_.second);
        }

        bool operator <(const TopPairing &other ) const{
            return indexes_<other.indexes_;
        }

        static bool orderByScore(const TopPairing & a, const TopPairing &  b){
            return a.score_ > b.score_;
        }

        friend std::ostream& operator<< (std::ostream &out, const TopPairing & topPairing){
            out<<topPairing.score_<<" "<<topPairing.indexes_.first<<" "<<topPairing.indexes_.second;

            return out;
        }

        std::pair<Genetics::Id::Snp, Genetics::Id::Snp> indexes_;
        float score_;

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
            inline TopSnpList(Genetics::Id::Snp numberSnps) :
                topK_(numberSnps),
                cutoff_(0.0),
                currentPartners_(std::vector<std::pair<Genetics::Id::Snp, float>>(numberSnps, { -1,0.0f })),
                testCounter_({0,0}),
                insertSincePrefix_(0)
            {
                pairwiseSignificanceCounts_.fill(0);
            }

            bool attemptInsert(Genetics::Id::Snp snpIndex1, Genetics::Id::Snp snpIndex2, float score);

            //Gets the number of snp Indexes used
            Genetics::Id::Snp size() const
            {
                return currentPartners_.size();
            }

            float getCutoff() const
            {
                return cutoff_;
            }

            void incrementTestCounter(const TestCounter& tests);

            const TestCounter& getTestCounter() const
            {
                return testCounter_;
            }

            void calculateFormattedResults();

            const FormattedPairs& getPairs() const
            {
                return formattedPairs_;
            }

        private:
            void calculateReciprocalPairs();

            //Stores the most up to date best pairing for each SNP
            std::vector<std::pair<Genetics::Id::Snp, float>> currentPartners_;

            TestCounter testCounter_;

            float cutoff_;
            int topK_;

            /**
            * Holds the count of detected pairs with <index> significance
            */
            std::array<Genetics::Id::Snp, MAX_CUTOFF+1> pairwiseSignificanceCounts_;

            Genetics::Id::Snp  insertSincePrefix_;

            //Formatted results
            FormattedPairs formattedPairs_;

    };
}
