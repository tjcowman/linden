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
#include <tuple>

using ID_Snp = std::uint32_t; //Maxium number of input Snps
using ID_Genotype = std::uint8_t; //Maximum number of genotype values 0, 1, 2 ...
using ID_Sample = std::uint32_t; //Maximum number of input samples

namespace ID_Invalid{
	constexpr ID_Snp Snp = std::numeric_limits<ID_Snp>::max();
}


static const uint32_t ESTIMATED_LD_RANGE = 1000000;

struct Location {
    

    Location(uint32_t chromosome, uint32_t basePair) : chromosome_(chromosome), basePair_(basePair)
    {
       // chromosome_ = chromosome;
       // basePair_ = basePair;
    }

    //TODO: change to explicit value representing invalid, ex: numeric_limits<>::max()
    Location()
    {
        chromosome_ = static_cast<uint32_t>(-1);
        basePair_ = static_cast<uint32_t>(-1);
      }
/*
    Location(
        uint32_t chromosome = static_cast<uint32_t>(-1), 
        uint32_t basePair = static_cast<uint32_t>(-1)
    ) : 
        chromosome_(chromosome), 
        basePair_(basePair)
    { }
*/
    bool operator==(const Location& lhs)const {
        return std::tie(chromosome_, basePair_) == std::tie(lhs.chromosome_, lhs.basePair_);
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

    uint32_t chromosome_;
    uint32_t basePair_;
};


struct Locus {
    std::string id;
    Location location;

    static void to_serial(std::ostream& os, const Locus& e) {

        //Writes the id length and id chars
        size_t length = e.id.size();
      //  std::cerr << "TMP LENGTH WRIT " <<length<< std::endl;
        os.write(reinterpret_cast<const char*>(&length), sizeof(size_t));
        os.write(reinterpret_cast<const char*>(e.id.data()), length * sizeof(char) );

        //Writes the location struct 
        os.write(reinterpret_cast<const char*>(&e.location), sizeof(Location));
        //os.write(reinterpret_cast<const char*>(&length), sizeof(size_t));
    }

    //BUG: SIZE_T IS A HUGE ISSUE WAS NOT 8 BYTES ON MY PC
    static Locus from_serial(std::istream& is) {
        Locus e;
        size_t length;
        
        is.read(reinterpret_cast<char*>(&length), 8);//sizeof(size_t));

        e.id.resize(length);
        is.read(const_cast<char*>(e.id.data()), length*sizeof(char));

        is.read(reinterpret_cast<char*>(&e.location), sizeof(Location));

        return e;
    }

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
