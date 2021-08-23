#include "Graph.h"

#include <catch2/catch_test_macros.hpp>
#include <iostream>

using TestGraph = Graph<uint16_t, uint16_t>;






TEST_CASE() {
	TestGraph g1(0);
	TestGraph g2(1);
	TestGraph g4(4);
	
	g1.print(std::cout);
	std::cout << "\n\n";

	g2.print(std::cout);
	std::cout << "\n\n";

	TestGraph g3 = TestGraph::joinToRoot(3, g1, g2);
	g3.print(std::cout);
	std::cout << "\n\n";

	TestGraph g5 = TestGraph::joinToRoot(5, g4, g3);
	g5.print(std::cout);
	std::cout << "\n";

	REQUIRE(1 == 1);
}
