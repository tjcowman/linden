/**
 * @author Tyler Cowman
 *
 * Acts as a wrapper for a collection of LD-Trees and associated functionality
 * such as merging and testing the set of trees.
 */


#pragma once

#include "Snp.hpp"
#include "SnpSet.h"
#include "LDTree.h"
#include "Locus.hpp"
#include "TopSnpList.h"

#include <limits.h>
#include <vector>
#include <fstream>

namespace Linden::Core
{
    //The number of ldtrees each tree will compare against for merging
    static const Genetics::Id::Snp MERGE_SEARCH_DISTANCE = 10;

    class LDForest
    {
    public:
        //TODO implement to not require a num snps by correctly sizing the topSnplist based on internal data (NOTE NEEEDS THE MAX POSSIBLE INDEX NOT FILTERED NUMBER)
        LDForest();
        LDForest( Genetics::Id::Snp numSnps);
        LDForest( SnpSet& snpSet, Genetics::Id::Snp numSnps); //TODO: Make only require SnpSet

        size_t size() const;
        bool operator==(const LDForest& lhs) const;

        void mergeTrees(double maxUnkownFraction);
        void testTrees(int maxThreadUsage);
        //void writeResults(const std::vector<Locus>& infoMatrix, Args& args);
        void writeResults(const std::vector<Linden::Genetics::Locus>& infoMatrix, const std::string& output);

        static void to_serial(std::ostream& os, const LDForest& e);
        static LDForest from_serial(std::istream& is);

    private:
        size_t mergeTreeIteration(float unknownFraction);

        std::vector<LDTree> ldtrees_;
        TopSnpList topSnpList_;
    };
}