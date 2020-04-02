#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include "Mesh.h"
#include <queue>
#include <functional>
#include <utility>
//#include <chrono>

using namespace std;

struct Functor
{
	bool operator()(const pair<float, int> &x, const pair<float, int> &y)
    {
		// First entries in the pairs are distances.
		return x.first < y.first;
	}
};

class Dijkstra
{
public:
	Dijkstra() {};
	vector<int> GetShortestPath(Mesh *mesh, int sourceIdx, int destIdx);
	pair<vector<float>, vector<int>> DijkstraMinHeap(const Mesh *mesh, int sourceIdx);
	
private:
	void SortVector(vector<pair<float,int>> &v);
	void GetPath(const vector<int> prev, int v, vector<int> &result);
};


#endif