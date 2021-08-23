/**
 * @author Tyler Cowman
 *
 * Class representing a single LD-Tree. 
 * New implementation in order to utlize non-full binary trees of SNPs
 * The SNP nodes are seperated from the genome locations to reduce the memory footprint when calculating contingecy tables.
 */

#pragma once

#include "Snp.h"
#include "TopSnpList.h"
#include "Graph.h"

#include <iostream>
#include <vector>


class LDForest;

class LDTree
{
    friend LDForest;
public:
    LDTree(const Snp& snp, const Location& location);
 
    LDTree( LDTree& t1,  LDTree& t2);

    bool empty()const;
    size_t size()const;

    const Snp& getRoot()const;
    const Snp& getLeft(ID_Snp id)const;
    const Snp& getRight(ID_Snp id)const;

    ID_Sample computeDifferences(const LDTree& other)const;
    bool validMerge(const LDTree& other, ID_Sample maxDiff)const;

    void clear();

    void epistasisTest(const LDTree& other)const;



private:

    ID_Snp root_; //Root node as the graph is being utilized as a directed tree
    Graph<Snp, ID_Snp> snps_;
    std::vector<Location> locations_; //Stores the genomic location of snps in this tree  //TODO: reimplement this 


    TopSnpList* topSnpList_;

    

   // std::vector<Snp> nodes_;
   // std::vector<Location> genomeLocations_;
};
