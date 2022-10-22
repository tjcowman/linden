
#include <array>
#include <charconv>
#include <string>
#include <map>

#include "InputParser.hpp"
#include "Location.hpp"

std::ifstream openFileChecked(std::string filepath){
    std::ifstream file(filepath);

    if (!file.is_open()){
        std::cerr << "file <" << filepath << "> not opened\n";
    }
    return file;

}

//TODO: Improve input validation ex: char at end of numeric values
std::vector<Linden::Genetics::Locus> parseLoci(std::ifstream ifs){

    if (!ifs.is_open())
    {
        return std::vector<Linden::Genetics::Locus>();
    }

    std::vector<Linden::Genetics::Locus> ret;
    std::map<std::string, uint32_t> chromosomeEncoder;

    //determine delimiter between tab and space by peeking the first line
    std::string lineBuffer;

    // Get current position
    auto len = ifs.tellg();
    std::getline(ifs, lineBuffer);
    std::string delimiter = lineBuffer.find("\t") != std::string::npos ? "\t": " " ;
    // Return to position before "Read line".
    ifs.seekg(len, std::ios_base::beg);

    while (std::getline(ifs, lineBuffer)) {
        std::array<uint32_t, 4> delims{0,0,0,0};
        for (uint8_t i = 0; i < 3; ++i) {
            auto nextDelim = lineBuffer.find(delimiter, delims[i]+1);
            delims[i+1] = (nextDelim != std::string::npos) ? nextDelim : lineBuffer.size();
        }

        uint32_t chromosome;
        uint32_t location;

        auto res = std::from_chars(lineBuffer.data()+delims[1]+1, lineBuffer.data() + delims[2], chromosome,10);
        if (res.ec != std::errc{}) { //must be represented by some non integral type
            std::cerr << "chromsome read error\n";
            return std::vector<Linden::Genetics::Locus>();
            
            std::string chr = lineBuffer.substr(delims[1] + 1, delims[2] - delims[1]);
            if (chromosomeEncoder.count(chr) == 0) {
                chromosomeEncoder[chr] = chromosomeEncoder.size()+23; //human numeric chromosomes from 1->22
            }
            chromosome = chromosomeEncoder[chr];
        }

        res = std::from_chars(lineBuffer.data() + delims[2]+1, lineBuffer.data() + delims[3], location);
        if (res.ec != std::errc{}) {
            std::cerr << "location read error\n";
            return std::vector<Linden::Genetics::Locus>();
        }
        ret.emplace_back(Linden::Genetics::Locus{ lineBuffer.substr(delims[0], delims[1] - delims[0]), Linden::Genetics::Location(chromosome, location) });
    }
    ifs.close();

    return ret;
}

GenotypeMatrix parseGenotypes(std::ifstream ifs) {

    if (!ifs.is_open())
    {
        return GenotypeMatrix();
    }

    GenotypeMatrix ret;

    std::string lineBuffer;
    ID_Sample width = 0;

    while (std::getline(ifs, lineBuffer)) {
        if (width == 0) {
            width = lineBuffer.size();
        }else if(width != lineBuffer.size()){
            std::cerr << "input matrix not rectangular" << "\n";
            return GenotypeMatrix();
        }

        for (const auto& e : lineBuffer) {
            ID_Genotype val = static_cast<ID_Genotype>(e - '0');
            if (val <= 2) {
                ret.data.push_back(val);
            }else{
                std::cerr << "input matrix has invalid genotype" << "\n";
                return GenotypeMatrix();
            }
        }
    }
    ifs.close();

    ret.width = width;
    ret.height = ret.data.size() / width;
    return ret;
}

