
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <omp.h>
#include <random>
#include <vector>

#include "argParser.h"
#include "InputParser.hpp"
#include "LDForest.h"
#include "Snp.hpp"
#include "SnpSet.h"
#include "TopSnpList.h"

struct Args
{
    std::string loci;
    std::string controls;
    std::string cases;

    //if a serialized snpset provided, dont need the other files
    std::filesystem::path snpSet = ""; 

    std::string output = "";

    int permuteSamples = 0;

    int maxThreads = 1;
    double maxUnknown = .1;
    double minMAF = 0.0;
    size_t maxMS = 1;

    void printSummaryRelevant(std::ofstream& ofs)
    {
        ofs << maxUnknown << "\t"
            << minMAF << "\t" << maxMS << "\t"
            << permuteSamples;
    }
};

namespace Linden
{ 
    // obtain a time-based seed:
    unsigned  randSeed = std::chrono::system_clock::now().time_since_epoch().count();
    struct Log
    {
        Genetics::Id::Snp snps_;
        Genetics::Id::Sample controls_;
        Genetics::Id::Sample cases_;

        //Varies based on trial
        Genetics::Id::Snp mafRemoved_;
        Genetics::Id::Snp marginalSignificanceRemoved_;

        Genetics::Id::Snp passingSnps_;
        Genetics::Id::Snp mergedTreesFormed_;

        void printDimensions(std::ofstream& ofs)
        {
            ofs << snps_ << "\t" << cases_ << "\t" << controls_;
        }

        void printSubsetDimensions(std::ofstream& ofs)
        {
            ofs << mafRemoved_ << "\t" << marginalSignificanceRemoved_ << "\t" << mergedTreesFormed_;
        }
    };

    void testData(Args& args){
        std::vector<Genetics::Locus> loci;
        Genetics::GenotypeMatrix cases;
        Genetics::GenotypeMatrix controls;

        Core::SnpSet snps;
        if (args.snpSet.empty())
        {
            //The controls and cases expect 0,1,2 chars
            #pragma omp parallel sections num_threads( std::min(3, args.maxThreads ) )
            {
                #pragma omp section
                loci = parseLoci(openFileChecked(args.loci));
                #pragma omp section
                cases = parseGenotypes(openFileChecked(args.cases));
                #pragma omp section
                controls = parseGenotypes(openFileChecked(args.controls));
            }
            std::clog << "files read" << "\n";

            //Set the Snp static size values
            //TODO REPLACE WITH SETTING FROM SNPSET AND CHANGE API IN SNP TO TAKE A SNPDIMENSION
            Core::Snp::setDimensions(controls.width, cases.width);

            snps = Core::SnpSet(loci, controls, cases);
        }
        /*else
        {
            std::cout<<"Reading from "<<args.snpSet<<std::endl;

            if(!std::filesystem::exists(args.snpSet))
            {
                std::cerr<<args.snpSet<<" does not exist"<<std::endl;
                return;
            }

            std::ifstream is(args.snpSet, std::ios::binary);

            if(!is.is_open())
            {
                std::cerr<<"failed to open file"<<std::endl;
                return;
            }

            snps = Core::SnpSet::from_serial(is);

            Core::Snp::setDimensions(snps.getDimensions().numControls_, snps.getDimensions().numCases_);

            std::cout<<"File Read "<<args.snpSet<<std::endl;
        }*/

        //Record the initial size of the dataset read in
        //TODO Rework this log class and try to remove
        Log log;
        log.snps_ = snps.size();//loci.size();
        log.cases_ = snps.getDimensions().numCases_;// cases.width;
        log.controls_ = snps.getDimensions().numControls_;// controls.width;


        //perform filtering logic on the input snps based on single locus measures
        std::clog << "filtering SNPs" << "\n";
        std::clog << "\tinitial: " << log.snps_ << "\n";

        log.mafRemoved_ = snps.remove_if([args](const Core::Snp& e){return e.computeMinorAlleleFrequency() < args.minMAF; });
        std::clog << "\tremoved minor allele frequency: " << log.mafRemoved_ << "\n";

        log.marginalSignificanceRemoved_ = snps.remove_if([args](const Core::Snp& e) {return e.marginalTest() > Statistics::GetChi2Table()[args.maxMS]; });
        std::clog << "\tremoved marginal significance: " << log.marginalSignificanceRemoved_ << "\n";

        //Initialize the ldforest, note that currently the loci size needs to refer to the range of possible indexes not how many post filterd SNPs there are
        //This is due to the implementation of TopSnpList
        Core::LDForest ldforest(snps, snps.getSizeUnfiltered() );

        if(ldforest.size() > 1)
        {
            ldforest.mergeTrees(args.maxUnknown);
            log.mergedTreesFormed_ = ldforest.size();
            ldforest.testTrees(args.maxThreads);
            ldforest.writeResults(snps.getLoci(), args.output);
        }
    }    
}

int main(int argc, char *argv[])
{
    srand(time(NULL));

    //Bring up help menu when no args, help, or h are first argument
    if( (argc == 1) || ((std::string)argv[1] == "--help") || ((std::string)argv[1] == "-h") )
    {
        std::cout<<"linden v 1.0"<<std::endl;
        return 0;
    }

    ARGLOOP(,
        ARG(maxThreads, stol)
        ARG(loci,)
        ARG(controls,)
        ARG(cases,)
        ARG(snpSet,)
        ARG(output,)
        ARG(permuteSamples, stoi)
        ARG(maxUnknown, stof)
        ARG(minMAF, stof)
        ARG(maxMS, stoi)
    )

    Linden::testData(args);

    return 0;
}
