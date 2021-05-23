#include "CommonStructs.h"


#include <string>
#include <charconv>

std::ifstream openFileChecked(std::string filepath){
    std::ifstream file(filepath);

    if (!file.is_open()){
        std::cerr << "file " << filepath << " not opened\n";
        exit(1);
    }
    else{
        return file;
    }
}

//TODO: Improve input validation ex: char at end of numeric values
std::vector<Locus> parseLoci(std::ifstream& ifs){
    std::vector<Locus> ret;
    const std::string delimiter = "\t";

    std::string lineBuffer;
    while (std::getline(ifs, lineBuffer)) {
   
        std::array<uint32_t, 4> delims{0};
        for (uint8_t i = 0; i < 3; ++i) {
            delims[i+1] = lineBuffer.find(delimiter, delims[i]+1);
        }

        uint32_t chromosome;
        uint32_t location;

        auto res = std::from_chars(lineBuffer.data()+delims[1]+1, lineBuffer.data() + delims[2], chromosome);
        if (res.ec != std::errc{}) {
            std::cerr << "chromsome read error\n";
            exit(1);
        }
        res = std::from_chars(lineBuffer.data() + delims[2]+1, lineBuffer.data() + delims[3], location);

        if (res.ec != std::errc{}) {
            std::cerr << "location read error\n";
            exit(1);
        }

        ret.emplace_back(Locus{ lineBuffer.substr(delims[0], delims[1]-delims[0]), chromosome, location });
    }

    return ret;
}

GenotypeMatrix parseGenotypes(std::ifstream& ifs) {
    GenotypeMatrix ret;

    std::string lineBuffer;
    uint32_t width = 0;

    while (std::getline(ifs, lineBuffer)) {
        if (width == 0) {
            width = lineBuffer.size();
        }else if(width != lineBuffer.size()){
            std::cerr << "input matrix not rectangular" << "\n";
            exit(1);
        }

        for (const auto& e : lineBuffer) {
            uint8_t val = e - '0';
            if (val <= 2) {
                ret.data.push_back(val);
            }else{
                std::cerr << "input matrix has invalid genotype" << "\n";
                exit(1);
            }
        }
    }
    ifs.close();

    ret.width = width;
    ret.height = ret.data.size() / width;
    return ret;
}