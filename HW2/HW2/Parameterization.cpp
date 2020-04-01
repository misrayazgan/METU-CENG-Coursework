#include "Parameterization.h"

Parameterization::Parameterization(WeightEnum w)
{
	weight = w;
}

bool Parameterization::isMouthVertex(Mesh *mesh, int v)
{
	float *coords = mesh->verts[v]->coords;
	if(coords[0] >= -19 && coords[0] <= 27 && 
		coords[1] >= -37 && coords[1] <= -33 &&
		coords[2] >= 4 && coords[2] <= 24)
	{
		return true;
	}
	return false;
}

set<int> Parameterization::FindBoundaryVertices(Mesh *mesh)
{
	set<int> boundaryVertices;
	set<int> outerBoundaryVertices;
	vector<FreeEdge *> boundaryEdges;
	vector<FreeEdge *> outerBoundaryEdges;
	isBoundary.resize(mesh->verts.size(), false);

	// First Dim: smaller valued vertex
	// Second Dim: larger valued vertex
	vector<vector<int>> edgeCount(mesh->verts.size(), vector<int>(mesh->verts.size(), 0));

	for(int k = 0; k < mesh->tris.size(); k++)
	{
		Triangle *tri = mesh->tris[k];
		int i, j;

		// 1 and 2
		i = min(tri->v1i, tri->v2i);
		j = max(tri->v1i, tri->v2i);
		edgeCount[i][j]++;

		// 2 and 3
		i = min(tri->v2i, tri->v3i);
		j = max(tri->v2i, tri->v3i);
		edgeCount[i][j]++;

		// 1 and 3
		i = min(tri->v1i, tri->v3i);
		j = max(tri->v1i, tri->v3i);
		edgeCount[i][j]++;
	}

	// Boundary edges are the edges which belong to only one face.
	for(int i = 0; i < edgeCount.size(); i++)
	{
		for(int j = 0; j < edgeCount[i].size(); j++)
		{
			if(edgeCount[i][j] == 1)
			{
				if(isMouthVertex(mesh, i) == false && isMouthVertex(mesh, j) == false)
				{
					boundaryVertices.insert(i);
					boundaryVertices.insert(j);

					FreeEdge *edge = new FreeEdge(i, j, mesh->computeEdgeLength(i, j));
					boundaryEdges.push_back(edge);
				}
			}
		}
	}

	// Store inner(smaller) and outer(larger) boundary
	set<int> boundaryVerts1;
	set<int> boundaryVerts2 = boundaryVertices;
	// Find first boundary edge 
	FreeEdge *firstEdge = boundaryEdges[0];
	int startV = firstEdge->v1i;
	int prevV = firstEdge->v1i;
	int currentV = firstEdge->v2i;
	boundaryVerts2.erase(startV);
	boundaryVerts2.erase(currentV);
	// Find next v
	int v = FindSuccessiveVertex(boundaryEdges, prevV, currentV);
	boundaryVerts1.insert(firstEdge->v1i);
	boundaryVerts1.insert(firstEdge->v2i);
	while(startV != v)
	{
		boundaryVerts1.insert(v);
		boundaryVerts2.erase(v);
		prevV = currentV;
		currentV = v;
		v = FindSuccessiveVertex(boundaryEdges, prevV, currentV);
	}

	// The outer boundary is the one with higher number of vertices.
	if(boundaryVerts1.size() > boundaryVerts2.size())
	{
		outerBoundaryVertices = boundaryVerts1;
	}
	else
	{
		outerBoundaryVertices = boundaryVerts2;
	}

	// Fill isBoundary to use later.
	set<int>::iterator it = outerBoundaryVertices.begin();
	for(; it != outerBoundaryVertices.end(); ++it)
	{
		int v = *it;
		isBoundary[v] = true;
	}

	return outerBoundaryVertices;
}

// Returns the vertex coming after currentV, uses prevV to compare.
int Parameterization::FindSuccessiveVertex(vector<FreeEdge *> edges, int prevV, int currentV)
{
	for(int i = 0; i < edges.size(); i++)
	{
		int v1i = edges[i]->v1i;
		int v2i = edges[i]->v2i;

		if(v1i == currentV && v2i != prevV)
		{
			return v2i;
		}
		if(v2i == currentV && v1i != prevV)
		{
			return v1i;
		}
	}
}

float Parameterization::FindAngle(float *v1coords, float *v2coords)
{
	float dotPro12 = v1coords[0] * v2coords[0] + v1coords[1] * v2coords[1] + v1coords[2] * v2coords[2];
	float dotPro11 = v1coords[0] * v1coords[0] + v1coords[1] * v1coords[1] + v1coords[2] * v1coords[2];
	float dotPro22 = v2coords[0] * v2coords[0] + v2coords[1] * v2coords[1] + v2coords[2] * v2coords[2];

	float result = acos(dotPro12 / sqrt(dotPro11 * dotPro22));
	return result;
}

float* Parameterization::GetVector(float* src, float* dst)
{
	float *vec = new float[3];
	vec[0] = dst[0] - src[0];
	vec[1] = dst[1] - src[1];
	vec[2] = dst[2] - src[2];

	return vec;
}

vector<FreeVertex *> Parameterization::ParamUniform(Mesh *mesh, vector<FreeVertex *> circleVertices)
{
	int numberOfVertices = mesh->verts.size();
	Eigen::MatrixXf W = Eigen::MatrixXf::Zero(numberOfVertices, numberOfVertices);
	Eigen::VectorXf bx(numberOfVertices);
	Eigen::VectorXf by(numberOfVertices);
	vector<FreeVertex *> resultVertices;

	// Fill W.
	for(int i = 0; i < numberOfVertices; i++)
	{
		for(int j = 0; j < numberOfVertices; j++)
		{
			// On the diagonal and boundary vertex
			if(i == j && isBoundary[i] == true)
			{
				W(i, j) = 1;
			}
			else if(i != j && mesh->isNeighbor(i, j) && isBoundary[i] == false)
			{
				W(i, j) = 1;
			}
		}
	}

	// On the diagonal and non-boundary vertex
	for(int i = 0; i < numberOfVertices; i++)
	{
		if(isBoundary[i] == false)
		{
			// Total value of the row
			W(i, i) = -1 * W.rowwise().sum()(i);
		}
	}

	// Fill bx and by.
	for(int i = 0; i < numberOfVertices; i++)
	{
		float *coords = new float[3];
		coords[0] = mesh->verts[i]->coords[0];
		coords[1] = mesh->verts[i]->coords[1];
		coords[2] = 0;

		if(isBoundary[i] == true)
		{
			// For each boundary vertex, find the closest circle vertex.
			int closestIdx = FindClosestCircleVertex(coords, circleVertices);
			FreeVertex *closestCircleV = circleVertices[closestIdx];
			//circleVertices.erase(circleVertices.begin() + closestIdx);
			bx(i) = closestCircleV->coords[0];
			by(i) = closestCircleV->coords[1];
		}
		else
		{
			// For the non-boundary vertices bx and by are zero.
			bx(i) = 0;
			by(i) = 0;
		}
	}

	// xx and xy store the coordinates of the new vertices which lie on the circle.
	// xx and xy are ordered with respect to the global idx of the vertices.
	Eigen::VectorXf xx = W.lu().solve(bx);
	Eigen::VectorXf xy = W.lu().solve(by);


	for(int i = 0; i < numberOfVertices; i++)
	{
		float *coords = new float[3];
		coords[0] = xx(i);
		coords[1] = xy(i);
		coords[2] = 0;
		FreeVertex *v = new FreeVertex(coords);
		resultVertices.push_back(v);
	}

	return resultVertices;
}

vector<FreeVertex *> Parameterization::ParamHarmonic(Mesh *mesh, vector<FreeVertex *> circleVertices)
{
	int numberOfVertices = mesh->verts.size();
	Eigen::MatrixXf W = Eigen::MatrixXf::Zero(numberOfVertices, numberOfVertices);
	Eigen::VectorXf bx(numberOfVertices);
	Eigen::VectorXf by(numberOfVertices);
	vector<FreeVertex *> resultVertices;

	// Fill W.
	// For each triangle, consider each edge.
	for(int i = 0; i < mesh->tris.size(); i++)
	{
		Triangle *tri = mesh->tris[i];
		float *v1coords = mesh->verts[tri->v1i]->coords;
		float *v2coords = mesh->verts[tri->v2i]->coords;
		float *v3coords = mesh->verts[tri->v3i]->coords;

		// 1 and 2 :: Vectors from 3 to 1 and 3 to 2
		float *vector1 = GetVector(v3coords, v1coords);
		float *vector2 = GetVector(v3coords, v2coords);
		float angle = FindAngle(vector1, vector2);

		if(isBoundary[tri->v1i] == false)
		{
			W(tri->v1i, tri->v2i) += (1.0 / tan(angle)) / 2;
		}
		if(isBoundary[tri->v2i] == false)
		{
			W(tri->v2i, tri->v1i) += (1.0 / tan(angle)) / 2;
		}		

		// 2 and 3 :: Vectors from 1 to 2 and 1 to 3
		vector1 = GetVector(v1coords, v2coords);
		vector2 = GetVector(v1coords, v3coords);
		angle = FindAngle(vector1, vector2);

		if(isBoundary[tri->v2i] == false)
		{
			W(tri->v2i, tri->v3i) += (1.0 / tan(angle)) / 2;
		}
		if(isBoundary[tri->v3i] == false)
		{
			W(tri->v3i, tri->v2i) += (1.0 / tan(angle)) / 2;
		}		

		// 1 and 3 :: Vectors from 2 to 1 and 2 to 3
		vector1 = GetVector(v2coords, v1coords);
		vector2 = GetVector(v2coords, v3coords);
		angle = FindAngle(vector1, vector2);

		if(isBoundary[tri->v1i] == false)
		{
			W(tri->v1i, tri->v3i) += (1.0 / tan(angle)) / 2;
		}
		if(isBoundary[tri->v3i] == false)
		{
			W(tri->v3i, tri->v1i) += (1.0 / tan(angle)) / 2;
		}
	}
	
	for(int i = 0; i < numberOfVertices; i++)
	{
		if(isBoundary[i] == true)
		{
			W(i, i) = 1;
		}
	}

	for(int i = 0; i < numberOfVertices; i++)
	{
		if(isBoundary[i] == false)
		{
			// Total value of the row
			W(i, i) = -1 * W.rowwise().sum()(i);
		}
	}

	// Fill bx and by.
	for(int i = 0; i < numberOfVertices; i++)
	{
		float *coords = new float[3];
		coords[0] = mesh->verts[i]->coords[0];
		coords[1] = mesh->verts[i]->coords[1];
		coords[2] = 0;

		if(isBoundary[i] == true)
		{
			// For each boundary vertex, find the closest circle vertex.
			int closestIdx = FindClosestCircleVertex(coords, circleVertices);
			FreeVertex *closestCircleV = circleVertices[closestIdx];
			bx(i) = closestCircleV->coords[0];
			by(i) = closestCircleV->coords[1];
		}
		else
		{
			// For the non-boundary vertices bx and by are zero.
			bx(i) = 0;
			by(i) = 0;
		}
	}

	// xx and xy store the coordinates of the new vertices which lie on the circle.
	// xx and xy are ordered with respect to the global idx of the vertices.
	Eigen::VectorXf xx = W.lu().solve(bx);
	Eigen::VectorXf xy = W.lu().solve(by);

	for(int i = 0; i < numberOfVertices; i++)
	{
		float *coords = new float[3];
		coords[0] = xx(i);
		coords[1] = xy(i);
		coords[2] = 0;
		FreeVertex *v = new FreeVertex(coords);
		resultVertices.push_back(v);
	}

	//writeMatrix(mesh, W);

	return resultVertices;
}


vector<FreeVertex *> Parameterization::ParamMeanValue(Mesh *mesh, vector<FreeVertex *> circleVertices)
{
	int numberOfVertices = mesh->verts.size();
	Eigen::MatrixXf W = Eigen::MatrixXf::Zero(numberOfVertices, numberOfVertices);
	Eigen::VectorXf bx(numberOfVertices);
	Eigen::VectorXf by(numberOfVertices);
	vector<FreeVertex *> resultVertices;

	// Fill W.
	// For each triangle, consider each edge.
	for(int i = 0; i < mesh->tris.size(); i++)
	{
		Triangle *tri = mesh->tris[i];
		float *v1coords = mesh->verts[tri->v1i]->coords;
		float *v2coords = mesh->verts[tri->v2i]->coords;
		float *v3coords = mesh->verts[tri->v3i]->coords;

		// 1 and 2 :: Vectors from 1 to 3 and 1 to 2
		float *vector1 = GetVector(v1coords, v3coords);
		float *vector2 = GetVector(v1coords, v2coords);
		float angle = FindAngle(vector1, vector2);
		float dist = FindDistance(v1coords, v2coords);

		if(isBoundary[tri->v1i] == false)
		{
			W(tri->v1i, tri->v2i) += (tan(angle / 2) / (2 * dist));
		}
		if(isBoundary[tri->v2i] == false)
		{
			W(tri->v2i, tri->v1i) += (tan(angle / 2) / (2 * dist));
		}		

		// 2 and 3 :: Vectors from 2 to 1 and 2 to 3
		vector1 = GetVector(v1coords, v2coords);
		vector2 = GetVector(v1coords, v3coords);
		angle = FindAngle(vector1, vector2);
		dist = FindDistance(v2coords, v3coords);

		if(isBoundary[tri->v2i] == false)
		{
			W(tri->v2i, tri->v3i) += (tan(angle / 2) / (2 * dist));
		}
		if(isBoundary[tri->v3i] == false)
		{
			W(tri->v3i, tri->v2i) += (tan(angle / 2) / (2 * dist));
		}		

		// 1 and 3 :: Vectors from 1 to 2 and 1 to 3
		vector1 = GetVector(v2coords, v1coords);
		vector2 = GetVector(v2coords, v3coords);
		angle = FindAngle(vector1, vector2);
		dist = FindDistance(v1coords, v3coords);

		if(isBoundary[tri->v1i] == false)
		{
			W(tri->v1i, tri->v3i) += (tan(angle / 2) / (2 * dist));
		}
		if(isBoundary[tri->v3i] == false)
		{
			W(tri->v3i, tri->v1i) += (tan(angle / 2) / (2 * dist));
		}
	}
	
	for(int i = 0; i < numberOfVertices; i++)
	{
		if(isBoundary[i] == true)
		{
			W(i, i) = 1;
		}
	}

	for(int i = 0; i < numberOfVertices; i++)
	{
		if(isBoundary[i] == false)
		{
			// Total value of the row
			W(i, i) = -1 * W.rowwise().sum()(i);
		}
	}

	// Fill bx and by.
	for(int i = 0; i < numberOfVertices; i++)
	{
		float *coords = new float[3];
		coords[0] = mesh->verts[i]->coords[0];
		coords[1] = mesh->verts[i]->coords[1];
		coords[2] = 0;

		if(isBoundary[i] == true)
		{
			// For each boundary vertex, find the closest circle vertex.
			int closestIdx = FindClosestCircleVertex(coords, circleVertices);
			FreeVertex *closestCircleV = circleVertices[closestIdx];
			bx(i) = closestCircleV->coords[0];
			by(i) = closestCircleV->coords[1];
		}
		else
		{
			// For the non-boundary vertices bx and by are zero.
			bx(i) = 0;
			by(i) = 0;
		}
	}

	// xx and xy store the coordinates of the new vertices which lie on the circle.
	// xx and xy are ordered with respect to the global idx of the vertices.
	Eigen::VectorXf xx = W.lu().solve(bx);
	Eigen::VectorXf xy = W.lu().solve(by);

	for(int i = 0; i < numberOfVertices; i++)
	{
		float *coords = new float[3];
		coords[0] = xx(i);
		coords[1] = xy(i);
		coords[2] = 0;
		FreeVertex *v = new FreeVertex(coords);
		resultVertices.push_back(v);
	}

	return resultVertices;
}

void Parameterization::writeMatrix(Mesh *mesh, Eigen::MatrixXf W)
{ 
	FILE* file = fopen("miso_output.txt", "w");

	if(file != NULL)
	{
		for(int i = 0; i < mesh->verts.size(); i++)
		{
			for(int j = 0; j < mesh->verts.size(); j++)
			{
				fprintf(file, "%.6g ", W(i, j));
			}
			fprintf(file, "\n");
		}

		fclose(file);
	}
} 

vector<FreeVertex *> Parameterization::GetCircleVertices(int numberOfCircleV)
{
	vector<FreeVertex *> circleVertices;

	// Center is (0,0,0)
	float *centerCoords = new float[3];
	centerCoords[0] = 0;
	centerCoords[0] = 0;
	centerCoords[0] = 0;
	FreeVertex *center = new FreeVertex(centerCoords);

	float radius = 10;
	float arcAngle = 360.0 / numberOfCircleV;
	
	for(int i = 0; i < numberOfCircleV; i++)
	{
		float angle = (arcAngle * i) * PI / 180;

		float *coords = new float[3];
		coords[0] = radius * sin(angle);
		coords[1] = radius * cos(angle);
		coords[2] = 0;

		FreeVertex *v = new FreeVertex(coords);
		circleVertices.push_back(v);
	}

	return circleVertices;
}

// Find the closest circle point to the source vertex
int Parameterization::FindClosestCircleVertex(float *srcCoords, vector<FreeVertex *> circleVertices)
{
	float inf = numeric_limits<float>::infinity();
	float minDistance = inf;
	FreeVertex *closestCircleV;
	int closestIdx = 0;

	for(int i = 0; i < circleVertices.size(); i++)
	{
		float dist = FindDistance(srcCoords, circleVertices[i]->coords);
		if(dist < minDistance)
		{
			minDistance = dist;
			closestCircleV = circleVertices[i];
			closestIdx = i;
		}
	}

	return closestIdx;
}

float Parameterization::FindDistance(float* coords1, float* coords2)
{
	float dist = sqrt(pow((coords2[0] - coords1[0]), 2) + pow((coords2[1] - coords1[1]), 2) + pow((coords2[2] - coords1[2]), 2));
	return dist;
}

vector<FreeVertex *> Parameterization::DiskParameterization(Mesh *mesh)
{
	vector<FreeVertex *> resultVertices;

	// Get boundary vertices
	set<int> boundaryVertices = FindBoundaryVertices(mesh);

	// Map boundary vertices to disk
	vector<FreeVertex *> circleVertices = GetCircleVertices(boundaryVertices.size());

	// Map non-boundary vertices
	if(weight == UNIFORM)
	{
		resultVertices = ParamUniform(mesh, circleVertices);
	}
	else if(weight == HARMONIC)
	{
		resultVertices = ParamHarmonic(mesh, circleVertices);
	}
	else if(weight == MEAN_VALUE)
	{
		resultVertices = ParamMeanValue(mesh, circleVertices);
	}

	return resultVertices;
}

//vector<FreeVertex *> Parameterization::ClosedMeshParam(Mesh *mesh)
//{
//}
//
//vector<FreeVertex *> Parameterization::SphereParam(Mesh *mesh, vector<int> boundaryVertices)
//{
//
//}