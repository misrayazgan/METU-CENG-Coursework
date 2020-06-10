#pragma once

#include "Mesh.h"

class ConvexHull
{
private:
	vector<Edge *> hullEdges;
public:
	void FindConvexHull(Mesh *mesh);
};


