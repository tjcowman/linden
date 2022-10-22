
#pragma once

#include "Types.hpp"

namespace Linden::Genetics
{
    struct Location
    {
        static const uint32_t ESTIMATED_LD_RANGE = 1000000;

        inline Location(uint32_t chromosome, uint32_t basePair) :
            chromosome_(chromosome),
            basePair_(basePair)
        { }

        //TODO: change to explicit value representing invalid, ex: numeric_limits<>::max()
        inline Location() :
            chromosome_(static_cast<uint32_t>(-1)),
            basePair_(static_cast<uint32_t>(-1))
        { }

        inline bool operator==(const Location& lhs) const
        {
            return std::tie(chromosome_, basePair_) == std::tie(lhs.chromosome_, lhs.basePair_);
        }

        inline friend std::ostream& operator<<(std::ostream& os, const Location& e)
        {
            os << e.chromosome_ << "\t" << e.basePair_;
            return os;
        }

        inline bool inLinkageDisequilibrium(const Location& other) const
        {
            if (chromosome_ != other.chromosome_)
            {
                return false;
            }
            //gets the absolute difference between two unsigned integral numbers
            ID_Snp dist = basePair_ > other.basePair_ ? basePair_ - other.basePair_ : other.basePair_ - basePair_;
            return(dist <= ESTIMATED_LD_RANGE);
        }

        uint32_t chromosome_;
        uint32_t basePair_;
    };
}