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
    int maxMS=1;
    
    
    void printSummaryRelevant(std::ofstream & ofs)
    {
        ofs <<maxUnknown<<"\t"
            <<minMAF<<"\t"<<maxMS<<"\t"
            <<permuteSamples;
    }
    
};

struct Locus {
    std::string id;
    uint32_t chromosome;
    uint32_t location;

   friend std::ostream& operator<<(std::ostream& os, const Locus& e) {
        os << e.id << "\t" << e.chromosome << "\t" << e.location;
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
    //uint32_t height() { return data.size() / width; }

};

struct DatasetSizeInfo
{
    long snps_;
    long cases_;
    long controls_;
    
    //Varies based on trial
    long mafRemoved_;
    long marginalSignificanceRemoved_;
    
    long passingSnps_;
    long mergedTreesFormed_;
    
    void printDimensions(std::ofstream & ofs)
    {
        ofs<<snps_<<"\t"<<cases_<<"\t"<<controls_;
    }
    
    void printSubsetDimensions(std::ofstream & ofs)
    {
        ofs<<mafRemoved_<<"\t"<<marginalSignificanceRemoved_<<"\t"<<mergedTreesFormed_;
    }
    
};


#endif //COMMON_STRUCTS_H
