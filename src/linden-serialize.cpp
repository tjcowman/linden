

#include <iostream>
#include <vector>
#include <omp.h>

#include "Types.h"
#include "argParser.h"
#include "InputParser.h"
#include "CommonStructs.h"
#include "Snp.h"
#include "SnpSet.h"
#include "LDForest.h"

struct Args {
    std::string loci;
    std::string controls;
    std::string cases;

    std::string output = "";

    //std::string 

    int maxThreads;
};


int main(int argc, char* argv[]) {

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



    std::vector<Locus> loci;
    GenotypeMatrix cases;
    GenotypeMatrix controls;


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

    SnpSet snps(loci, controls, cases);
    //LDForest ldforest(snps, loci.size());

    if (args.output != "") {
        std::ofstream ofs(args.output, std::ofstream::binary);
        SnpSet::to_serial(ofs, snps);
        //LDForest::to_serial(ofs, ldforest);
        ofs.close();
    }
    else { //print to standard out
       // LDForest::to_serial(std::cout, ldforest);
        SnpSet::to_serial(std::cout, snps);
    }

}