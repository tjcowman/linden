

#include <iostream>
#include <vector>
#include <omp.h>

#include "argParser.h"
#include "InputParser.hpp"
#include "Snp.hpp"
#include "SnpSet.h"
#include "LDForest.h"

struct Args {
    std::string loci;
    std::string controls;
    std::string cases;

    std::string output = "";
    int maxThreads;
};


int main(int argc, char* argv[])
{

    //Bring up help menu when no args, help, or h are first argument
    if ((argc == 1) || ((std::string)argv[1] == "--help") || ((std::string)argv[1] == "-h")) {
        std::cout << "linden v 1.0" << std::endl;
        return 0;
    }


    ARGLOOP(,
        ARG(maxThreads, stol)
        ARG(loci, )
        ARG(controls, )
        ARG(cases, )
        ARG(output, )
    )



    std::vector<Linden::Genetics::Locus> loci;
    Linden::Genetics::GenotypeMatrix cases;
    Linden::Genetics::GenotypeMatrix controls;


    //The controls and cases expect 0,1,2 chars
    #pragma omp parallel sections num_threads( std::min(3, args.maxThreads ) )
    {
        #pragma omp section
        loci = Linden::parseLoci(Linden::openFileChecked(args.loci));
        #pragma omp section
        cases = Linden::parseGenotypes(Linden::openFileChecked(args.cases));
        #pragma omp section
        controls = Linden::parseGenotypes(Linden::openFileChecked(args.controls));
    }

    Linden::Core::SnpSet snps(loci, controls, cases);
    //LDForest ldforest(snps, loci.size());

    if (args.output != "") {
        std::ofstream ofs(args.output, std::ofstream::binary);
        Linden::Core::SnpSet::to_serial(ofs, snps);
        //LDForest::to_serial(ofs, ldforest);
        ofs.close();
    }
    else { //print to standard out
       // LDForest::to_serial(std::cout, ldforest);
        Linden::Core::SnpSet::to_serial(std::cout, snps);
    }

}