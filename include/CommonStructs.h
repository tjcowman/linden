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



struct Args{
    std::string fileInfo;
    std::string fileControls;
    std::string fileCases;
    
    std::string output="out";
    
    int permuteSamples=0;
    
    int maxThreads=1;
    float maxUnknown=.1;
    float minMAF=0.0;
    int maxMS=1;
    
    
    void printSummaryRelevant(std::ofstream & ofs)
    {
        ofs <<maxUnknown<<"\t"
            <<minMAF<<"\t"<<maxMS<<"\t"
            <<permuteSamples;
    }
    
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
