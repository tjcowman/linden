/**
 * @author Tyler Cowman
 * 
 * Class representing a single LD-Tree. The SNP nodes are seperated from the genome locations 
 * to reduce the memory footprint when calculating contingecy tables.
 */

#ifndef LDTREE_H
#define LDTREE_H

#include "Snp.h"
#include "TopSnpList.h"

#include <iostream>
#include <vector>

#define ESTIMATED_LD_RANGE 1000000

/*struct GenomeLocation
{
    GenomeLocation(char chromosome, int basePair)
    {
        chromosome_ = chromosome;
        basePair_ = basePair;
    }
    
    GenomeLocation(const GenomeLocation & l1, const GenomeLocation & l2)
    {
        chromosome_ = '!';
        basePair_ = -1;
    }
    
    bool inLinkageDisequilibrium(const GenomeLocation & other)const
    {
        if(chromosome_ != other.chromosome_)
            return false;
        else if(abs(basePair_ - other.basePair_) > ESTIMATED_LD_RANGE)
            return false;
        else
            return true;
            
    }
    
    char chromosome_;
    int basePair_;
};
*/

class LDTree
{
    public:
      
      //  LDTree(const Snp & snp, char chromosome, int basePair);
        LDTree(const Snp& snp, const Location& location);

        LDTree(const LDTree & cpy);
        
        LDTree(const LDTree & t1, const LDTree & t2);
    
        bool empty()const;
        size_t size()const;
        
        const Snp & getRoot()const;
        
        int computeDifferences(const LDTree & other)const;
        
        void clear();
        
        void epistasisTest(const LDTree & other, TopSnpList & topSnpList)const;
        void epistasisTestNoTrees(const LDTree & other, TopSnpList & topSnpList)const;
        
        friend std::ostream& operator<< (std::ostream &out, const LDTree & ldtree);
        
    private:
        
        std::vector<Snp> nodes_;
        //std::vector<GenomeLocation> genomeLocations_;
        std::vector<Location> genomeLocations_;

};
#endif //LDTREE_H
