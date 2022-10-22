
#pragma once

#include <iostream>
#include <map>
#include <vector>

#include "Graph.h"
#include "Location.hpp"
#include "Snp.hpp"
#include "TopSnpList.h"

////////////////////////////////////////////////////////////////////////////////
//! @class LDTree
//!
//! Class representing a single LD-Tree. New implementation in order to utlize 
//! non-full binary trees of SNPs. The SNP nodes are seperated from the genome 
//! locations to reduce the memory footprint when calculating contingecy tables.
////////////////////////////////////////////////////////////////////////////////
class LDTree
{
public:
    // Constructor
    inline LDTree(const Snp& snp, const Linden::Genetics::Location& location, TopSnpList* topSnpList = nullptr) :
        snps_(Graph<Snp, ID_Snp>(snp)),
        locations_(Graph<Linden::Genetics::Location, ID_Snp>(location)),
        topSnpList_( topSnpList),
        root_(0)
    { }

    ////////////////////////////////////////////////////////////////////////////
    //! Constructs a new LDTree by merging the root nodes of two other
    //! LDTrees.
    //!
    //! @param t1 The left tree node.
    //! @param t2 The right tree node.
    ////////////////////////////////////////////////////////////////////////////
    inline LDTree( LDTree& t1,  LDTree& t2) :
        root_(0),
        locations_(Graph<Linden::Genetics::Location, ID_Snp>::joinToRoot(
            Linden::Genetics::Location(),
            t1.locations_,
            t2.locations_
        )),
        snps_(Graph<Snp, ID_Snp>::joinToRoot(
            Snp(t1.getRoot(), t2.getRoot()),
            t1.snps_,
            t2.snps_
        )),
        topSnpList_(nullptr)
    {
        if(t1.topSnpList_ == t2.topSnpList_)
        {
            topSnpList_ = t1.topSnpList_;
        }
        else
        {
            std::cerr<<"subTrees point to different output lists"<<std::endl;
        }
    }


    const Snp& getRoot() const;

    inline std::vector<ID_Snp> getChildren(ID_Snp i ) const 
    {
        return snps_.getOutgoingIds(i);
    }

    inline auto getChildren2(ID_Snp i) const
    {
        return snps_.getOutgoingIdsItStart(i);
    }

    bool isLeaf(ID_Snp i) const;


    ID_Sample computeDifferences(const LDTree& other) const;
    bool validMerge(const LDTree& other, ID_Sample maxDiff) const;


    void epistasisTest(const LDTree& other) const;
    
    // Clears the data contained in the tree.
    inline void clear() 
    {
        root_ = ID_Invalid::Snp;
        snps_.clear();
        locations_.clear();
        topSnpList_ = nullptr;
    }

    // Gets if the tree is empty.
    inline bool empty() const
    {
        return snps_.empty();
    }

    // Gets the number of snp nodes in the tree.
    inline size_t size() const
    {
        return snps_.size();
    }

    ////////////////////////////////////////////////////////////////////////
    //! Gets the snp with the given Id.
    //!
    //! @param id The Snp id to lookup.
    //! @returns The Snp.
    /////////////////////////////////////////////////////////////////////////
    inline const Snp& getSnp(ID_Snp id) const
    {
        return snps_.getElement(id);
    }

    // Gets the output list used by this tree.
    const TopSnpList* getTopSnpList()
    {
        return topSnpList_;
    }

    // Equivalence operator compares all but the output list.
    inline bool operator==(const LDTree& lhs) const 
    {
        return (std::tie(root_, snps_, locations_) == 
                std::tie(lhs.root_, lhs.snps_, lhs.locations_));
    }

    static void to_serial(std::ostream& os, const LDTree& e);
    static LDTree from_serial(std::istream& is);

private:

    //! Root node of tree represented by the graph structure.
    ID_Snp root_;
    //! Genotype data for snps in the tree.
    Graph<Snp, ID_Snp> snps_;
    //! Location of snps in the tree.
    Graph<Linden::Genetics::Location, ID_Snp> locations_; 
    //! List used to track snp pairs, shared between all trees.
    TopSnpList* topSnpList_;
};
