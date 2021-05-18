#include <iostream>
#include <omp.h>
#include <vector>
#include <fstream>
#include <map>
#include <time.h>

#include "argParser.h"
#include "CommonStructs.h"

#include "Snp.h"
#include "TopSnpList.h"
#include "LDForest.h"


#include "Matrix.h"
#include "MatrixMath.h"
#include "TsvParser.h"


//Used for determining the chi2 correseponding to a pvalue of 1*10^-i for single locus significance
//Input as an integer denoting the -log10 signficance up to 6
const static std::array<float, 7> chi2DegreesFreedomTable{0, 2.71, 6.63, 10.82, 15.14, 19.57, 24.87};


void testData(Args& args)
{
    TsvParser T;
    
    //The infoFile expects plaintext 
    T.setDelimiter('\t');

    
    Matrix<char> casesMatrix; 
    Matrix<char> controlsMatrix;
    Matrix<std::string> lociMatrix = T.parse(args.loci);

    //The controls and cases expect 0,1,2 chars
    T.setDelimiter(' ');
    #pragma omp parallel sections num_threads( std::min(2, args.maxThreads ) )
    {        
        #pragma omp section
        casesMatrix = T.parseByteWise(args.cases);
        #pragma omp section
        controlsMatrix = T.parseByteWise(args.controls);    
    }
    
    //If any of the files were read incorrectly exit
    if(controlsMatrix.size() == 0 || casesMatrix.size() == 0 || lociMatrix.size() == 0)
    {
        std::cerr<<"input file error"<<"\n";
        exit(1);
    }
    
    
    //Record the initial size of the dataset read in
    DatasetSizeInfo datasetSizeInfo;
    datasetSizeInfo.snps_ = lociMatrix.dim(0);
    datasetSizeInfo.cases_ = casesMatrix.dim(1);
    datasetSizeInfo.controls_ = controlsMatrix.dim(1);

    datasetSizeInfo.mafRemoved_ = 0;
    datasetSizeInfo.marginalSignificanceRemoved_ = 0;

    //Call snp constructors to create bitwise snp representations
    std::vector<Snp> snps;
    

    if( args.permuteSamples != 1) 
    {
        for(long i= 0; i<lociMatrix.dim(0); ++i)
        {
            snps.push_back(Snp(i, controlsMatrix.getRow(i), casesMatrix.getRow(i)));
        }
    }
    else
    {
        std::pair<Matrix<char>, Matrix<char> > permutedMatrixes = MatrixMath::permuteColumns(controlsMatrix, casesMatrix);

        for(long i= 0; i<lociMatrix.dim(0); ++i)
            snps.push_back(Snp(i, permutedMatrixes.first.getRow(i), permutedMatrixes.second.getRow(i)));
    }

    LDForest ldforest( lociMatrix.dim(0) , controlsMatrix.dim(1), casesMatrix.dim(1), lociMatrix.dim(0));


    //Only create trees from snps with a high enough MAF and low enough marginal significance
    for(int i=0; i<snps.size(); ++i)
    {  
        if(snps[i].computeMinorAlleleFrequency() <  args.minMAF )
        {
            datasetSizeInfo.mafRemoved_++; 
        }
        else if(snps[i].marginalTest() > chi2DegreesFreedomTable[args.maxMS])
        {
            datasetSizeInfo.marginalSignificanceRemoved_++;
        }
        else
        {
            if(isdigit(lociMatrix.a(snps[i].getIndex(),1)[0]));
            {
                ldforest.insert(snps[i], (char)stoi(lociMatrix.a(snps[i].getIndex(),1)), stoi(lociMatrix.a(snps[i].getIndex(),2))); 
            }   
        }    
    }

    datasetSizeInfo.passingSnps_ = ldforest.size();
    
    std::clog<<"filtering SNPs"<<"\n";
    std::clog<<"\tinitial: "<<datasetSizeInfo.snps_<<"\n";
    std::clog<<"\tremoved marginal significance: "<<datasetSizeInfo.marginalSignificanceRemoved_<<"\n";
    std::clog<<"\tremoved minor allele frequency: "<<datasetSizeInfo.mafRemoved_<<"\n";
    std::clog<<"\tremaining: "<<ldforest.size()<<"\n";
    
    
    
    if(ldforest.size() > 1)
    {    
        ldforest.mergeTrees( args.maxUnknown, datasetSizeInfo);
        datasetSizeInfo.mergedTreesFormed_ = ldforest.size();
        ldforest.testTrees(args.maxThreads);
        ldforest.writeResults(lociMatrix, args, datasetSizeInfo);
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
        ARG(output,)
        ARG(permuteSamples, stoi)
        ARG(maxUnknown, stof)
        ARG(minMAF, stof)
        ARG(maxMS, stoi)
    );
 
    
  
    testData(args);
}
