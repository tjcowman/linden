/**
 * @author Tyler Cowman
 *  
 * Acts as a wrapper for a collection of LD-Trees and associated functionality
 * such as merging and testing the set of trees.
 */


#pragma once

#include "CommonStructs.h"
#include "Snp.h"
#include "SnpSet.h"

//Include the original Tree version for testing
#ifdef oldTree
#include "LDTree.h"
#else
#include "LDTree2.h"
#endif

#include "TopSnpList.h"
#include <limits.h>
#include <vector>
#include <fstream>

//The number of ldtrees each tree will compare against for merging
static const ID_Snp MERGE_SEARCH_DISTANCE = 10;

class LDForest{
    public:

        //TODO implement to not require a num snps by correctly sizing the topSnplist based on internal data (NOTE NEEEDS THE MAX POSSIBLE INDEX NOT FILTERED NUMBER)
        LDForest();
        LDForest( ID_Snp numSnps);
        LDForest( SnpSet& snpSet, ID_Snp numSnps); //TODO: Make only require SnpSet
        
        size_t size()const;
        bool operator==(const LDForest& lhs)const;

        void mergeTrees(double maxUnkownFraction);
        void testTrees(int maxThreadUsage);   
        //void writeResults(const std::vector<Locus>& infoMatrix, Args& args);
        void writeResults(const std::vector<Locus>& infoMatrix, const std::string& output);

        static void to_serial(std::ostream& os, const LDForest& e);
        static LDForest from_serial(std::istream& is);

    private:
        size_t mergeTreeIteration(float unknownFraction);
        
        std::vector<LDTree> ldtrees_;
        TopSnpList topSnpList_;

           

};
