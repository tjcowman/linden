#include <iostream>
#include <omp.h>
#include <vector>
#include <fstream>
#include <map>
#include <time.h>
#include <chrono>
#include <random>
#include <algorithm>

#include "argParser.h"
#include "CommonStructs.h"

#include "Snp.h"
#include "SnpSet.h"
#include "TopSnpList.h"
#include "LDForest.h"

#include "InputParser.h"




//Used for determining the chi2 correseponding to a pvalue of 1*10^-i for single locus significance
//Input as an integer denoting the -log10 signficance up to 6
const static std::array<float, 7> chi2DegreesFreedomTable{0.0f, 2.71f, 6.63f, 10.82f, 15.14f, 19.57f, 24.87f};

// obtain a time-based seed:
unsigned  randSeed = std::chrono::system_clock::now().time_since_epoch().count();

void testData(Args& args){
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
    std::clog << "files read" << "\n";

    

    //Call snp constructors to create bitwise snp representations
    std::vector<Snp> snps;
    if (args.permuteSamples != 1) {
        for (size_t i = 0; i < loci.size(); ++i) {
            snps.push_back(Snp(i, controls, cases));
        }
    }
    else {

        std::vector<uint32_t> indexOrder;
        indexOrder.reserve(cases.width + controls.width);
        for (ID_Sample i = 0; i < cases.width + controls.width; ++i) {
            indexOrder.push_back(i);
        }

        std::shuffle(indexOrder.begin(), indexOrder.end(), std::default_random_engine(randSeed));
        /*
        for (int i = 0; i < ordering.size(); ++i) {
            if (i < m1.dim(1)){
                if (ordering[i] < m1.dim(1)){
                    retVal.first.addColumn(m1.getColumn(ordering[i]));
                }
                else{
                    retVal.first.addColumn(m2.getColumn(ordering[i] - m1.dim(1)));
                }

            }
            else {
                if (ordering[i] < m1.dim(1)){
                    retVal.second.addColumn(m1.getColumn(ordering[i]));
                }
                else{
                    retVal.second.addColumn(m2.getColumn(ordering[i] - m1.dim(1)));
                }
            }
        }

        // std::pair<Matrix<char>, Matrix<char> > permutedMatrixes = MatrixMath::permuteColumns(controlsMatrix, casesMatrix);

        //  for (long i = 0; i < lociMatrix.dim(0); ++i)
        //    snps.push_back(Snp(i, permutedMatrixes.first.getRow(i), permutedMatrixes.second.getRow(i)));
        */
    }
    
    //Record the initial size of the dataset read in
    //TODO Rework this log class and try to remove
    Log log{};
    log.snps_ = loci.size();
    log.cases_ = cases.width;
    log.controls_ = controls.width;

    //Set the Snp static size values
    //TODO REPLACE WITH SETTING FROM SNPSET AND CHANGE API IN SNP TO TAKE A SNPDIMENSION
    Snp::setDimensions(controls.width, cases.width);
   

    SnpSet snps_N(loci, controls, cases);

    //perform filtering logic on the input snps based on single locus measures
    std::clog << "filtering SNPs" << "\n";
    std::clog << "\tinitial: " << log.snps_ << "\n";

    log.mafRemoved_ = snps_N.remove_if([args](const Snp& e){return e.computeMinorAlleleFrequency() < args.minMAF; });
    std::clog << "\tremoved minor allele frequency: " << log.mafRemoved_ << "\n";

    log.marginalSignificanceRemoved_ = snps_N.remove_if([args](const Snp& e) {return e.marginalTest() > chi2DegreesFreedomTable[args.maxMS]; });
    std::clog << "\tremoved marginal significance: " << log.marginalSignificanceRemoved_ << "\n";


    //Initialize the ldforest, note that currently the loci size needs to refer to the range of possible indexes not how many post filterd SNPs there are
    //This is due to the implementation of TopSnpList
    //BUG PROBABLY I THINK THIS LOCATION INDEX IS WRONG WHEN FILTERED
    LDForest ldforest(snps_N, loci.size() );


    
    
    if(ldforest.size() > 1){    
        ldforest.mergeTrees( args.maxUnknown);
        log.mergedTreesFormed_ = ldforest.size();

        //LDForest::to_serial(std::cout, ldforest);
        //exit(0);

        ldforest.testTrees(args.maxThreads);
        ldforest.writeResults(loci, args);
    }

    


}

int main(int argc, char *argv[]){
    srand(time(NULL));
    
    //Bring up help menu when no args, help, or h are first argument
    if( (argc == 1) || ((std::string)argv[1] == "--help") || ((std::string)argv[1] == "-h") ) {
        std::cout<<"linden v 1.0"<<std::endl;
        return 0;
    }
    
    
    ARGLOOP(,
        ARG(maxThreads, stol)
        ARG(loci,)
        ARG(controls,)
        ARG(cases,)
        ARG(output,)
        ARG(permuteSamples, stoi)
        ARG(maxUnknown, stof)
        ARG(minMAF, stof)
        ARG(maxMS, stoi)
    )
 
    
  
    testData(args);
}
