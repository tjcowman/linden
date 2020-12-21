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
#include "Matrix.h"


#include <limits.h>
#include <vector>
#include <fstream>


#define MERGE_SEARCH_DISTANCE 10

class LDForest
{
    public:
        LDForest(int topKSnps, int numberControlSamples, int numberCaseSamples, int numberSnps);
        
        void insert(const Snp & snp, char chromosome, int basePair);
        
        int size()const;
        
        void mergeTrees(float maxUnkownFraction, DatasetSizeInfo datasetSizeInfo);
        void testTrees(int maxThreadUsage);
        
        void writeResults(const Matrix<std::string>& infoMatrix, Args& args, DatasetSizeInfo datasetSizeInfo);
    
    private:
        int mergeTreeIteration(float unknownFraction, DatasetSizeInfo datasetSizeInfo);
        
        //Results output
        void writeGroundTruthList();
        
        std::vector<LDTree> ldtrees_;
        
        TopSnpList topSnpList_;
        
};
#endif //LDFOREST_H
