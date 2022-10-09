# Summary
LinDeN is a tool for computing statistically significantly interacting epistatic pairs of single nucleotide polymorphisms.

# Publication

https://doi.org/10.1093/nar/gkx505

# Installation

```bash
#Using the command line, navigate to the main linden directory and then create a build directory.
foo@bar:~$ mkdir build
#Enter the build directory.
foo@bar:~$ cd build
#Generate the build file using cmake.
foo@bar:~$ cmake ..
#Compile linden to create the binary in the build directory.
foo@bar:~$ make linden
#Or install linden, note this probably requires sudo access.
foo@bar:~$ make install
```

# Running
LinDen utilizes three files representing the set of SNPs, control genotypes and case genotypes as the primary input. Note that these three files are required and are currently the only format readable by LinDeN. Examples are provided in the example directory.
```bash
foo@bar:~$ linden --loci <filepath> --controls <filepath> --cases <filepath>
```

## Loci
A three column tab seperated text file.


|SNP Name (string)|Chromosome (int)|Base Pair (int)|
|---|---|---|
SNP1	|1	|100000
SNP2	|1	|200000
...|...|...|


## Controls and Cases
A plaintext dense matrix of genotypes with rows corresponding to the SNPs in the loci file, and columns for the relevant number of samples. The controls and cases should be in seperate files.

```
00110100
00110101
...
```


# Other Arguments
```bash
#Filepath where the resulting cutoffPairs and reciprocalPairs will be created.
#If blank will print both standard out.
* --output

#Floating point value representing the maximum unknown genotype ratio in internal nodes of LD trees.
#Larger values are faster but result in a less accurate search.
#Roughly 0.45 is optimal, see the paper for details.
* --maxUnknown

#Filters out loci in which the minor allele frequency is too low for effective analysis.
#A commonly used threshold is 0.05
* --minMAF

#Takes an integer value from 1-6.
#Represents the maximum marginal significance of loci for consideration in pairwise testing as a -log10 p-value.
#For example, a value of 3 will filter out all loci with a marginal significance less than 0.001
* --maxMS

#Sets the maximum number of threads that will be utilized by LinDen.
#This should beset to the number of precessing units LinDen is intended to be run on.
#The runtime scaling with the numberof processing units is roughly linear.
* --maxThreads

#This can be set to 0 or 1 meaning permute and don’t permute respectively.
#In most cases this shoud be turned off.
#Permuting randomly assigns each sample (column) as a case or control.
#This is useful for permutation testing to obtain a measure of pairwise significances of loci in the input dataset under a null distribution.
* --permuteSamples

```
# Output
Creates a file for all pairs passing the final dynamic significance threshold. And a file for the smaller subset of reciprocally signficant pairs.

|chi2|	snp1|	snp2|	chr1|	chr2|	bp1|	bp2|
|---|---|---|---|---|---|--|
|3.2584|	SNP3|	SNP15|	1|	1|	300000|	1500000|
2.31468|	SNP5|	SNP16|	1|	1|	500000|	1600000|
