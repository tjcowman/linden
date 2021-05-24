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
#include "TopSnpList.h"
#include "LDForest.h"

#include "InputParser.h"




//Used for determining the chi2 correseponding to a pvalue of 1*10^-i for single locus significance
//Input as an integer denoting the -log10 signficance up to 6
const static std::array<float, 7> chi2DegreesFreedomTable{0, 2.71, 6.63, 10.82, 15.14, 19.57, 24.87};

// obtain a time-based seed:
unsigned randSeed = std::chrono::system_clock::now().time_since_epoch().count();

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
    Log log;
   // DatasetSizeInfo datasetSizeInfo;
    log.snps_ = loci.size();
    log.cases_ = cases.width;
    log.controls_ = controls.width;

    log.mafRemoved_ = 0;
    log.marginalSignificanceRemoved_ = 0;

    //LDForest ldforest(loci.size(), controls.width, cases.width, loci.size());
    LDForest ldforest(&log, loci.size());

    //Only create trees from snps with a high enough MAF and low enough marginal significance
    for(ID_Snp i=0; i<snps.size(); ++i){
       // std::cout << i << " " << snps[i].computeMinorAlleleFrequency() << std::endl;
        if(snps[i].computeMinorAlleleFrequency() <  args.minMAF ){
            log.mafRemoved_++;
        }
        else if(snps[i].marginalTest() > chi2DegreesFreedomTable[args.maxMS]){
            log.marginalSignificanceRemoved_++;
        }
        else{
          //  ldforest.insert(snps[i], loci[i].location.chromosome, loci[i].location.basePair);
            //LDTree l(snps[i], loci[i].chromosome, loci[i].location);
           // std::cout << loci[i] << "\n";
            ldforest.insert(LDTree(snps[i], loci[i].location));
           // std::cout << ldforest.size()<<"    ?\n";
        }    
    }

    log.passingSnps_ = ldforest.size();
    

    std::clog<<"filtering SNPs"<<"\n";
    std::clog<<"\tinitial: "<< log.snps_<<"\n";
    std::clog<<"\tremoved marginal significance: "<< log.marginalSignificanceRemoved_<<"\n";
    std::clog<<"\tremoved minor allele frequency: "<< log.mafRemoved_<<"\n";
    std::clog<<"\tremaining: "<<ldforest.size()<<"\n";
    
    
    
    if(ldforest.size() > 1){    
        ldforest.mergeTrees( args.maxUnknown);
        log.mergedTreesFormed_ = ldforest.size();
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
    );
 
    
  
    testData(args);
}
