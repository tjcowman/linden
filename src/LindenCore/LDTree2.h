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
#include <map>


class LDForest;

class LDTree
{
    friend LDForest;
public:
    LDTree();

    LDTree(const Snp& snp, const Location& location);
 
    LDTree( LDTree& t1,  LDTree& t2);

    bool operator==(const LDTree& lhs)const;

    bool empty()const;
    size_t size()const;

    const Snp& getSnp(ID_Snp i)const;
    const Snp& getRoot()const;

    std::vector<ID_Snp> getChildren(ID_Snp i)const;
    
    auto getChildren2(ID_Snp i) const
    {
        return snps_.getOutgoingIdsItStart(i);
    }

    bool isLeaf(ID_Snp i)const;


    ID_Sample computeDifferences(const LDTree& other)const;
    bool validMerge(const LDTree& other, ID_Sample maxDiff)const;

    void clear();

    void epistasisTest(const LDTree& other)const;

    static void to_serial(std::ostream& os, const LDTree& e);
    static LDTree from_serial(std::istream& is);

private:

    ID_Snp root_; //Root node as the graph is being utilized as a directed tree
    Graph<Snp, ID_Snp> snps_;
    Graph<Location, ID_Snp> locations_; //Stores the genomic location of snps in this tree  

    TopSnpList* topSnpList_;
};
