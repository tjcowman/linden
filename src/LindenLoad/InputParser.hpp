
#pragma once

#include <filesystem>
#include <fstream>
#include <vector>

#include "CommonStructs.h"

std::ifstream openFileChecked(std::string filepath);

//TODO: Improve input validation ex: char at end of numeric values
std::vector<Locus> parseLoci(std::ifstream ifs);

GenotypeMatrix parseGenotypes(std::ifstream ifs);
