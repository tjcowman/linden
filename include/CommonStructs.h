/**
 * @author Tyler Cowman
 * 
 * These structs serve as containers for the information about a given run of LinDen
 * and are primarily used to provide an organized location for this data when printing
 * to standard out or a file later.
 */


#ifndef COMMON_STRUCTS_H
#define COMMON_STRUCTS_H

#include <iostream>
#include <fstream>
#include <vector>

#include "Types.h"

struct Args{
    std::string loci;
    std::string controls;
    std::string cases;
    
    std::string output="";
    
    int permuteSamples=0;
    
    int maxThreads=1;
    double maxUnknown=.1;
    double minMAF=0.0;
    size_t maxMS=1;
     
    void printSummaryRelevant(std::ofstream & ofs)
    {
        ofs <<maxUnknown<<"\t"
            <<minMAF<<"\t"<<maxMS<<"\t"
            <<permuteSamples;
    }
};


//static const uint32_t ESTIMATED_LD_RANGE = 1000000;
static const uint32_t ESTIMATED_LD_RANGE = 1000000;
struct Location {
    uint32_t chromosome_;
    uint32_t basePair_;


    Location(uint32_t chromosome, uint32_t basePair)
    {
        chromosome_ = chromosome;
        basePair_ = basePair;
    }

    //TODO: change to explicit value representing invalid, ex: numeric_limits<>::max()
    Location()
    {
        chromosome_ = static_cast<uint32_t>(-1);
        basePair_ = static_cast<uint32_t>(-1);
      }

    friend std::ostream& operator<<(std::ostream& os, const Location& e) {
        os << e.chromosome_ << "\t" << e.basePair_;
        return os;
    }

    bool inLinkageDisequilibrium(const Location& other)const{
        if (chromosome_ != other.chromosome_) {
            return false;
        }
        //gets the absolute difference between two unsigned integral numbers
        ID_Snp dist = basePair_ > other.basePair_ ? basePair_ - other.basePair_ : other.basePair_ - basePair_;
        return(dist <= ESTIMATED_LD_RANGE);
    }
};

struct Locus {
    std::string id;
    Location location;

   friend std::ostream& operator<<(std::ostream& os, const Locus& e) {
        os << e.id << "\t" << e.location;
        return os;
    }
};

struct GenotypeMatrix {
    ID_Sample width;
    ID_Snp height;
    std::vector<ID_Genotype> data;

    std::vector<ID_Genotype>::const_iterator rowBegin(size_t row)const {
        return data.begin() + row * width;
    }

    std::vector<ID_Genotype>::const_iterator rowEnd(size_t row)const {
        return data.begin() + (row+1) * width;
    }
};

struct Log {
    ID_Snp snps_;
    ID_Sample controls_;
    ID_Sample cases_;
    
    //Varies based on trial
    ID_Snp mafRemoved_;
    ID_Snp marginalSignificanceRemoved_;

    ID_Snp passingSnps_;
    ID_Snp mergedTreesFormed_;

    void printDimensions(std::ofstream& ofs){
        ofs << snps_ << "\t" << cases_ << "\t" << controls_;
    }

    void printSubsetDimensions(std::ofstream& ofs){
        ofs << mafRemoved_ << "\t" << marginalSignificanceRemoved_ << "\t" << mergedTreesFormed_;
    }

};

#endif //COMMON_STRUCTS_H
