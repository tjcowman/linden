/**
 * @author Tyler Cowman
 *
 * Class representing a single LD-Tree. The SNP nodes are seperated from the genome locations
 * to reduce the memory footprint when calculating contingecy tables.
 */

#pragma once

#include "Snp.hpp"
#include "TopSnpList.h"

#include <iostream>
#include <vector>

#include "Location.hpp"

namespace Linden::Core
{

    class LDForest;

    class LDTree
    {
        friend LDForest;
    public:

        LDTree(const Snp& snp, const Genetics::Location& location, TopSnpList* topSnpList=nullptr) :
            nodes_(std::vector<Snp>{ snp }),
            genomeLocations_(std::vector<Genetics::Location>{ location }),
            topSnpList_(topSnpList)
        { }

        LDTree();

        // LDTree(const LDTree & cpy);
        LDTree( LDTree& t1,  LDTree& t2);

        bool operator==(const LDTree& lhs) const;

        bool empty() const;
        size_t size() const;

        const Snp& getRoot() const;

        Genetics::Id::Sample computeDifferences(const LDTree& other) const;
        bool validMerge(const LDTree& other, Genetics::Id::Sample maxDiff) const;

        void clear();

        void epistasisTest(const LDTree& other) const;


        static void to_serial(std::ostream& os, const LDTree& e);
        static LDTree from_serial(std::istream& is);

    private:

        TopSnpList* topSnpList_;

        std::vector<Snp> nodes_;
        std::vector<Genetics::Location> genomeLocations_;
    };
}
