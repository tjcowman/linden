/**
 * @author Tyler Cowman
 *  
 * Acts as a wrapper for a collection of LD-Trees and associated functionality
 * such as merging and testing the set of trees.
 */


#ifndef LDFOREST_H
#define LDFOREST_H

#include "CommonStructs.h"
#include "Snp.h"
#include "LDTree.h"
#include "TopSnpList.h"

#include <limits.h>
#include <vector>
#include <fstream>


static const ID_Snp MERGE_SEARCH_DISTANCE = 10;

class LDForest{
    public:
        LDForest(Log* log, ID_Snp numSnps);
        
        void insert( LDTree&& ldtree);
        
        size_t size()const;
        
        void mergeTrees(double maxUnkownFraction);
        void testTrees(int maxThreadUsage);   
        void writeResults(const std::vector<Locus>& infoMatrix, Args& args);
    
    private:
        size_t mergeTreeIteration(float unknownFraction);
        
        std::vector<LDTree> ldtrees_;
        TopSnpList topSnpList_;

        //Non-owning pointer to a log struct for storing statistics about the current run
        Log* log_;        

};
#endif //LDFOREST_H
