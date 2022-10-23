
#pragma once

#include <string>
#include <iostream>

#include "Location.hpp"
#include "Id.hpp"

namespace Linden::Genetics
{
    struct Locus
    {
        static void to_serial(std::ostream& os, const Locus& e)
        {
            //Writes the id length and id chars
            size_t length = e.id.size();
            os.write(reinterpret_cast<const char*>(&length), sizeof(size_t));
            os.write(reinterpret_cast<const char*>(e.id.data()), length * sizeof(char) );

            //Writes the location struct
            os.write(reinterpret_cast<const char*>(&e.location), sizeof(Location));
        }

        //BUG: SIZE_T IS A HUGE ISSUE WAS NOT 8 BYTES ON MY PC
        static Locus from_serial(std::istream& is)
        {
            Locus e;
            size_t length;

            is.read(reinterpret_cast<char*>(&length), 8);//sizeof(size_t));

            e.id.resize(length);
            is.read(const_cast<char*>(e.id.data()), length*sizeof(char));

            is.read(reinterpret_cast<char*>(&e.location), sizeof(Location));

            return e;
        }

        friend std::ostream& operator<<(std::ostream& os, const Locus& e)
        {
            os << e.id << "\t" << e.location;
            return os;
        }

        std::string id;
        Location location;
    };
}