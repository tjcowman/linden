#include "Types.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>

#pragma once




template<class T, class IT>
class Graph {
public:
	Graph();
	//Creates a new Graph of a single node containing a type of T and indexing type of IT
	Graph(T data);
	Graph(const Graph& cpy);

	bool operator==(const Graph& lhs)const;

	bool empty()const;
	size_t size()const;
	void clear();

	//Gets the element of type T from the passed index
	const T& getElement(IT id)const;

	std::vector<IT> getOutgoingIds(IT id)const;
	bool isTerminal(IT id)const;


	void print(std::ostream& os)const;


	//Increments the index for all vertices and edges, increasing the number of vertices
	//void addOffset(IT offset);

	//Creates a new graph representation by treating each graph as a tree rooted at index zero and joining them at the new root
	//Should probably move this out to the LDTree class as it is specific to trees and not graphs
	static Graph joinToRoot(T newRoot,  Graph& g1, Graph& g2 );

	static void to_serial(std::ostream& os, const Graph& e);
	static Graph from_serial(std::istream& is);

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
Graph<T, IT>::Graph(const Graph& cpy) : V(cpy.V), A(cpy.A), JA(cpy.JA), IA(cpy.IA){}

template<class T, class IT>
bool Graph<T, IT>::operator==(const Graph& lhs)const {
	return std::tie(V, A, JA, IA) == std::tie(lhs.V, lhs.A, lhs.JA, lhs.IA);

}

template<class T, class IT>
Graph<T, IT>::Graph(T data) {
	V = std::vector<T>{ data };
	A = std::vector<IT>();
	JA = std::vector<IT>();
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
std::vector<IT> Graph<T, IT>::getOutgoingIds(IT id)const{
	return std::vector<IT>(JA.begin() + IA[id], JA.begin() + IA[id + 1]);
}

template<class T, class IT>
bool Graph<T, IT>::isTerminal(IT id)const {
	return IA[id + 1] - IA[id] == 0;
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


template<class T, class IT>
void Graph<T, IT>::to_serial(std::ostream& os, const Graph& e) {

	//The T may be more complex than an integer indexing value
	IT num=e.V.size();
	os.write(reinterpret_cast<const char*>(&num), sizeof(IT));
	for (IT i = 0; i < num; ++i) {
		T::to_serial(os, e.V[i]);
	}

	vector_to_serial<IT, IT>(os, e.A);
	vector_to_serial<IT, IT>(os, e.JA);
	vector_to_serial<IT, IT>(os, e.IA);
}

template<class T, class IT>
Graph<T, IT>  Graph<T, IT>::from_serial(std::istream& is) {
	Graph e;

	IT num;
	is.read(reinterpret_cast<char*>(&num), sizeof(IT));
	e.V.reserve(num);
	for (IT i = 0; i < num; ++i) {
		e.V.push_back(T::from_serial(is));
	}
	e.A = vector_from_serial<IT, IT>(is);
	e.JA = vector_from_serial<IT,IT>(is);
	e.IA = vector_from_serial<IT, IT>(is);

	
	return e;
}


template<class T, class IT>
Graph<T,IT> Graph<T, IT>::joinToRoot(T newRoot, Graph& g1, Graph& g2) {
	Graph<T, IT> newGraph;
	IT offsetLeft = 1;
	IT offsetRight = offsetLeft + g1.V.size();

	//std::cout << "TS1 " << g1.V.size() << std::endl;
	//Move the node data
	newGraph.V = std::vector<T>{ newRoot };
	newGraph.V.insert(newGraph.V.end(), std::make_move_iterator(g1.V.begin()), std::make_move_iterator(g1.V.end()));
	newGraph.V.insert(newGraph.V.end(), std::make_move_iterator(g2.V.begin()), std::make_move_iterator(g2.V.end()));
	g1.V.clear();
	g2.V.clear();
	//std::cout <<"TS2 "<<g1.V.size() << std::endl;

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
	//auto offset = newGraph.IA.back();
	for ( auto it = g1.IA.begin() + 1; it != g1.IA.end(); ++it)
		//newGraph.IA.push_back(*it - (*(it - 1)) + offset);
		newGraph.IA.push_back(*it - (*(it - 1)) + newGraph.IA.back());

	//offset = newGraph.IA.back();
	for (auto it = g2.IA.begin() + 1; it != g2.IA.end(); ++it)
		newGraph.IA.push_back(*it - (*(it - 1)) + newGraph.IA.back());
		//newGraph.IA.push_back(*it - (*(it - 1)) + offset);

	//std::for_each(g1.IA.begin(), g1.IA.end(), [offsetLeft](auto& e) {e += offsetLeft; });
	//newGraph.IA.insert(newGraph.IA.end(), std::make_move_iterator(g1.IA.begin()), std::make_move_iterator(g1.IA.end()));
	//std::for_each(g2.IA.begin(), g2.IA.end(), [offsetRight](auto& e) {e += offsetRight; });
	//newGraph.IA.insert(newGraph.IA.end(), std::make_move_iterator(g2.IA.begin()), std::make_move_iterator(g2.IA.end()));

	return newGraph;
}