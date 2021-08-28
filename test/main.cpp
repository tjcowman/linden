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

	REQUIRE(g1.getOutgoingIds(0).size() == 0);
	REQUIRE(g1.isTerminal(0));

	g2.print(std::cout);
	std::cout << "\n\n";

	TestGraph g3 = TestGraph::joinToRoot(3, g1, g2);
	g3.print(std::cout);
	std::cout << "\n\n";

	auto ga(g3);
	auto gb(g3);
	auto gs1(g3);
	auto gs2(g3);
	auto gc = TestGraph::joinToRoot(0, ga, gb);
	gc.print(std::cout);
	std::cout<< "\n\n";

	TestGraph g5 = TestGraph::joinToRoot(5, g4, g3);
	g5.print(std::cout);
	std::cout << "\n\n";



	REQUIRE(g3.getOutgoingIds(0).size() == 2);
	REQUIRE(!g3.isTerminal(0));

	auto gb1(g5);
	auto gb2(g5);

	TestGraph gr = TestGraph::joinToRoot(7, gb1, gs1);
	gr.print(std::cout);
	std::cout << "\n\n";
	TestGraph gl = TestGraph::joinToRoot(7, gs2, gb2);
	gl.print(std::cout);
	std::cout << "\n";


	
}
