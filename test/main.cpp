#include "Graph.h"
#include "InputParser.h"
#include "Snp.h"
#include "argParser.h"
#include "CommonStructs.h"
#include "LDForest.h"

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <sstream>


//Integer type extension for allowing serialize and de-serialize with same api as Snps
struct IntType {
	IntType() {};
	IntType(uint16_t i): val(i){
		//std::cout << val << std::endl;
	}

	bool operator==(const IntType& lhs)const {
		return val == lhs.val;
	}

	friend std::ostream& operator<<(std::ostream& os, const IntType& e) {
		os << e.val;
		return os;
	}

	static IntType from_serial(std::istream& is) {
		IntType e;
		is.read(reinterpret_cast<char*>(&e.val), sizeof(IntType));
		return e;
	}

	static void to_serial(std::ostream& os, const IntType& e){
		os.write(reinterpret_cast<const char*>(&e.val), sizeof(IntType));
	}

	uint16_t val;
};


using TestGraph = Graph<IntType, uint16_t>;


//Test a serialization followed by a deserialization operation for equivalence on the passed object
template<class T>
bool testObject(const T& object) {
	std::stringstream res;
	T::to_serial(res, object);
	T objectNew = T::from_serial(res);

	return object == objectNew;
}

//Test for Graph structure serialization and deserialization
TEST_CASE() {
	TestGraph g1(0);
	TestGraph g2(1);
	REQUIRE(!(g1 == g2));

	TestGraph g3 = TestGraph::joinToRoot(3, g1, g2);

	//std::stringstream res;
	//TestGraph::to_serial(res, g3);
	//TestGraph g4 = TestGraph::from_serial(res);

	//REQUIRE(g3 == g4);
	REQUIRE(testObject(g3));

}

TEST_CASE() {
	
	std::vector<Locus> loci = parseLoci(openFileChecked("example/loci"));
	GenotypeMatrix cases = parseGenotypes(openFileChecked("example/cases"));
	GenotypeMatrix controls = parseGenotypes(openFileChecked("example/controls"));

	//Call snp constructors to create bitwise snp representations
	std::vector<Snp> snps;
	for (size_t i = 0; i < loci.size(); ++i) {
		snps.push_back(Snp(i, controls, cases));
	}
	Snp::setDimensions(controls.width, cases.width);

	REQUIRE(testObject(snps[0]));

	LDForest ldforest(nullptr, loci.size());

	for (ID_Snp i = 0; i < snps.size(); ++i) {
		ldforest.insert(LDTree(snps[i], loci[i].location));
	}

	//REQUIRE(testObject(ldforest));
	
}


TEST_CASE() {
	TestGraph g1(0);
	TestGraph g2(1);
	TestGraph g4(4);
	
//	g1.print(std::cout);
	//std::cout << "\n\n";

	REQUIRE(g1.getOutgoingIds(0).size() == 0);
	REQUIRE(g1.isTerminal(0));

	//g2.print(std::cout);
	//std::cout << "\n\n";

	TestGraph g3 = TestGraph::joinToRoot(3, g1, g2);
	//g3.print(std::cout);
	//std::cout << "\n\n";

	auto ga(g3);
	auto gb(g3);
	auto gs1(g3);
	auto gs2(g3);
	auto gc = TestGraph::joinToRoot(0, ga, gb);
	//gc.print(std::cout);
	//std::cout<< "\n\n";

	TestGraph g5 = TestGraph::joinToRoot(5, g4, g3);
	//g5.print(std::cout);
	//std::cout << "\n\n";



	REQUIRE(g3.getOutgoingIds(0).size() == 2);
	REQUIRE(!g3.isTerminal(0));

	auto gb1(g5);
	auto gb2(g5);

	TestGraph gr = TestGraph::joinToRoot(7, gb1, gs1);
	//gr.print(std::cout);
	//std::cout << "\n\n";
	TestGraph gl = TestGraph::joinToRoot(7, gs2, gb2);
//	gl.print(std::cout);
	//std::cout << "\n";


	
}
