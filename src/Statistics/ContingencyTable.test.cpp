
#include <chrono>
#include "catch2.hpp"

#include "ContingencyTable.hpp"

TEST_CASE("ContingencyTable")
{
    ContingencyTable<2> table;
    table.Data() = {
        34, 56, 90,
        22, 34, 56,
        56, 90
    };

    REQUIRE(abs(table.Chi2() - .0332) < .001);

}