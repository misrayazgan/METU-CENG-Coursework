#include "Mesh.h"
#include <queue>
#include <functional>
#include <utility>
//#include <chrono>

using namespace std;

enum DijkstraType
{
	DIJKSTRA_ARRAY = 1,
	DIJKSTRA_MINHEAP = 2,
	DIJKSTRA_FIBOHEAP = 3
};


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
	Dijkstra(DijkstraType type);
	vector<int> GetShortestPath(Mesh *mesh, int sourceIdx, int destIdx);
	vector<vector<float>> FindAllDistances(Mesh *mesh, bool write);
	pair<vector<float>, vector<int>> DijkstraArray(const Mesh *mesh, int sourceIdx);
	pair<vector<float>, vector<int>> DijkstraMinHeap(const Mesh *mesh, int sourceIdx);
	pair<vector<float>, vector<int>> DijkstraFiboHeap(const Mesh *mesh, int sourceIdx);
	
private:
	// Node structure for Fibonacci heap
	//struct Node {
	//	int numberOfChildren;
	//	Node *leftSib;       // Pointer to the left sibling
	//	Node *rightSib;      // Pointer to the right sibling
	//	Node *parent;
	//	Node *child;      // Pointer to the first child of a node
	//	bool marked;        // Whether the node is marked black or white
	//	float distance;   // Distance from the source node
	//	int idx;
	//	int degree;
	//};
	//Node *minPtr;
	//int numberOfFiboNodes;

	DijkstraType type;
	void SortVector(vector<pair<float,int>> &v);
	//bool Compare(const pair<float, int> &x, const pair<float, int> &y);
	void WriteMatrix(char* filename, vector<vector<float>> matrix);
	void GetPath(const vector<int> prev, int v, vector<int> &result);
};


