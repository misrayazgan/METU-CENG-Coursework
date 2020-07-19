#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>
#include <numeric>
#include <algorithm>
#include "Mesh.h"
#include "Dijkstra.h"


class Sampling
{
public:
	Sampling(int nSamples);
	//Vertex* FindFurthestVertex(vector<Vertex *> selectedVertices);
	vector<int> FPS(Mesh *mesh);
private:
	int numberOfSamples;
};


#endif