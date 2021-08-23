#include "Types.h"

#include <vector>
#include <iostream>
#include <algorithm>

template<class T, class IT>
class Graph {
public:
	Graph();
	//Creates a new Graph of a single node containing a type of T and indexing type of IT
	Graph(T data);

	bool empty()const;
	size_t size()const;
	void clear();

	//Gets the element of type T from the passed index
	const T& getElement(IT id)const;


	void print(std::ostream& os)const;


	//Increments the index for all vertices and edges, increasing the number of vertices
	//void addOffset(IT offset);

	//Creates a new graph representation by treating each graph as a tree rooted at index zero and joining them at the new root
	//Should probably move this out to the LDTree class as it is specific to trees and not graphs
	static Graph joinToRoot(T newRoot,  Graph& g1, Graph& g2 );

private:
	std::vector<T> V;
	std::vector<IT> A;
	std::vector<IT> JA;
	std::vector<IT> IA;
};

template<class T, class IT>
Graph<T, IT>::Graph() {

}

template<class T, class IT>
Graph<T, IT>::Graph(T data) {
	V = std::vector<T>{ data };
	//A = std::vector<IT>{ 1 }
	IA = std::vector<IT>{ 0,0 };
}

template<class T, class IT>
size_t Graph<T, IT>::size()const {
	return V.size();
}

template<class T, class IT>
const T& Graph<T, IT>::getElement(IT id)const{
	return V[id];
}

template<class T, class IT>
void Graph<T, IT>::clear() {
	V.clear();
	A.clear();
	JA.clear();
	IA.clear();
}

template<class T, class IT>
void Graph<T, IT>::print(std::ostream& os)const {
	os << "V: ";
	for (const auto& e : V) {
		os << e << " ";
	}
	os << "\n";

	os << "A: ";
	for (const auto& e : A) {
		os << e << " ";
	}
	os << "\n";

	os << "JA: ";
	for (const auto& e : JA) {
		os << e << " ";
	}
	os << "\n";

	os << "IA: ";
	for (const auto& e : IA) {
		os << e << " ";
	}
}


template<class T, class IT>
bool Graph<T, IT>::empty()const {
	return V.empty();
}


/*template<class T, class IT>
void Graph<T, IT>::addOffset(IT offset) {


}*/

template<class T, class IT>
Graph<T,IT> Graph<T, IT>::joinToRoot(T newRoot, Graph& g1, Graph& g2) {
	Graph<T, IT> newGraph;
	IT offsetLeft = 1;
	IT offsetRight = offsetLeft + g1.V.size();

	//Move the node data
	newGraph.V = std::vector<T>{ newRoot };
	newGraph.V.insert(newGraph.V.end(), std::make_move_iterator(g1.V.begin()), std::make_move_iterator(g1.V.end()));
	newGraph.V.insert(newGraph.V.end(), std::make_move_iterator(g2.V.begin()), std::make_move_iterator(g2.V.end()));

	//Expand the edge weights vector
	newGraph.A = std::vector<IT>{ 1,1 }; //Edge weights for the 2 new edges from root
	newGraph.A.insert(newGraph.A.end(), std::make_move_iterator(g1.A.begin()), std::make_move_iterator(g1.A.end()));
	newGraph.A.insert(newGraph.A.end(), std::make_move_iterator(g2.A.begin()), std::make_move_iterator(g2.A.end()));


	//Expand the column index vector
	//newGraph.JA = std::vector<IT>{ 0,0 };
	newGraph.JA = std::vector<IT>{ offsetLeft, offsetRight };
	std::for_each(g1.JA.begin(), g1.JA.end(), [offsetLeft](auto& e) {e += offsetLeft; });
	newGraph.JA.insert(newGraph.JA.end(), std::make_move_iterator(g1.JA.begin()), std::make_move_iterator(g1.JA.end()));
	std::for_each(g2.JA.begin(), g2.JA.end(), [offsetRight](auto& e) {e += offsetRight; });
	newGraph.JA.insert(newGraph.JA.end(), std::make_move_iterator(g2.JA.begin()), std::make_move_iterator(g2.JA.end()));

	//expand row indexes vector
	newGraph.IA = std::vector<IT>{ 0,2 };
	auto offset = newGraph.IA.back();
	for ( auto it = g1.IA.begin() + 1; it != g1.IA.end(); ++it)
		newGraph.IA.push_back(*it - (*(it - 1)) + offset);

	offset = newGraph.IA.back();
	for (auto it = g2.IA.begin() + 1; it != g2.IA.end(); ++it)
		newGraph.IA.push_back(*it - (*(it - 1)) + offset);

	//std::for_each(g1.IA.begin(), g1.IA.end(), [offsetLeft](auto& e) {e += offsetLeft; });
	//newGraph.IA.insert(newGraph.IA.end(), std::make_move_iterator(g1.IA.begin()), std::make_move_iterator(g1.IA.end()));
	//std::for_each(g2.IA.begin(), g2.IA.end(), [offsetRight](auto& e) {e += offsetRight; });
	//newGraph.IA.insert(newGraph.IA.end(), std::make_move_iterator(g2.IA.begin()), std::make_move_iterator(g2.IA.end()));

	return newGraph;
}