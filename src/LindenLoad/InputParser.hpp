
#pragma once

#include <filesystem>
#include <fstream>
#include <vector>


#include "Locus.hpp"
#include "GenotypeMatrix.hpp"

namespace Linden
{
    std::ifstream openFileChecked(std::string filepath);

    //TODO: Improve input validation ex: char at end of numeric values
    std::vector<Linden::Genetics::Locus> parseLoci(std::ifstream ifs);

    Genetics::GenotypeMatrix parseGenotypes(std::ifstream ifs);
}
