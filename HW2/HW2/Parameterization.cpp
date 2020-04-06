#include "Parameterization.h"

Parameterization::Parameterization(WeightEnum w, bool closed)
{
	weight = w;
	isClosed = closed;
}

set<int> Parameterization::FindBoundaryVertices(Mesh *mesh)
{
	set<int> boundaryVertices;
	set<int> outerBoundaryVertices;
	vector<FreeEdge *> boundaryEdges;
	vector<FreeEdge *> outerBoundaryEdges;

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
				boundaryVertices.insert(i);
				boundaryVertices.insert(j);

				FreeEdge *edge = new FreeEdge(i, j, mesh->computeEdgeLength(i, j));
				boundaryEdges.push_back(edge);
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

	SetBoundaryVertices(mesh, outerBoundaryVertices);

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

// If v is in cutVertices, return its index.
// Else return -1.
int Parameterization::GetIndex(const vector<int> cutVertices, int v)
{
	int result = find(cutVertices.begin(), cutVertices.end(), v) - cutVertices.begin();
	if(result < cutVertices.size())
	{
		return result; 
	}
	return -1;
}

void Parameterization::SetBoundaryVertices(Mesh *mesh, set<int> boundaryVertices)
{
	isBoundary.resize(mesh->verts.size(), false);

	// Fill isBoundary to use later.
	set<int>::iterator it = boundaryVertices.begin();
	for(; it != boundaryVertices.end(); ++it)
	{
		int v = *it;
		isBoundary[v] = true;
	}
}

vector<FreeVertex *> Parameterization::ParamUniform(Mesh *mesh, vector<FreeVertex *> circleVertices, vector<int> boundaryVerts)
{
	int numberOfVertices = mesh->verts.size();
	Eigen::SparseMatrix<float> W(numberOfVertices, numberOfVertices);
	Eigen::VectorXf bx(numberOfVertices);
	Eigen::VectorXf by(numberOfVertices);
	vector<FreeVertex *> resultVertices;

	// Fill W.
	for(int i = 0; i < numberOfVertices; i++)
	{
		// On the diagonal and boundary vertex
		if(isBoundary[i] == true)
		{
			W.insert(i, i) = 1;
		}
		else
		{
			// Not on the diagonal (i != j)
			// Total value of the row
			float sum = 0;
			vector<int> neighborVerts = mesh->verts[i]->vertList;
			for(int j = 0; j < neighborVerts.size(); j++)
			{
				W.insert(i, neighborVerts[j]) = 1;
				sum += 1;
			}
			// On the diagonal and non-boundary vertex
			W.insert(i, i) = -1 * sum;	// coeffref mi olacak ????
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
			FreeVertex *closestCircleV;
			if(isClosed == false)
			{
				int closestIdx = FindClosestCircleVertex(coords, circleVertices);
				closestCircleV = circleVertices[closestIdx];
				//circleVertices.erase(circleVertices.begin() + closestIdx);
			}
			else	// for closed mesh
			{
				// Find position of current v in boundaryVerts
				int idx = GetIndex(boundaryVerts, i);
				closestCircleV = new FreeVertex(circleVertices[idx]->coords);
			}
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
	Eigen::SparseLU<Eigen::SparseMatrix<float>, Eigen::AMDOrdering<int>> solver;
	solver.analyzePattern(W); 
	solver.factorize(W); 
	Eigen::VectorXf xx = solver.solve(bx);
	Eigen::VectorXf xy = solver.solve(by);

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
		else
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
		else
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

vector<int> Parameterization::FindPoleVertices(Mesh *mesh)
{
	vector<int> poles;
	int pole1 = 0;
	poles.push_back(pole1);

	Dijkstra *d = new Dijkstra();
	pair<vector<float>, vector<int>> result = d->DijkstraMinHeap(mesh, 0);
	vector<float> distances = result.first;
	int pole2 = max_element(distances.begin(), distances.end()) - distances.begin();
	poles.push_back(pole2);

	return poles;
}

void Parameterization::FindCutVertices(Mesh *mesh, vector<int> &cutVertices)
{
	// Find start and end points of the cut.
	vector<int> poles = FindPoleVertices(mesh);
	Dijkstra *d = new Dijkstra();
	cutVertices = d->GetShortestPath(mesh, poles[0], poles[1]);
}

// Given two triangles find the common edge vertices.
vector<int> Parameterization::FindCommonVertices(Mesh *mesh, int t1, int t2)
{
	Triangle *tri1 = mesh->tris[t1];
	Triangle *tri2 = mesh->tris[t2];
	vector<int> commonVerts;

	if(tri1->v1i == tri2->v1i || tri1->v1i == tri2->v2i || tri1->v1i == tri2->v3i)
	{
		commonVerts.push_back(tri1->v1i);
	}
	if(tri1->v2i == tri2->v1i || tri1->v2i == tri2->v2i || tri1->v2i == tri2->v3i)
	{
		commonVerts.push_back(tri1->v2i);
	}
	if(tri1->v3i == tri2->v1i || tri1->v3i == tri2->v2i || tri1->v3i == tri2->v3i)
	{
		commonVerts.push_back(tri1->v3i);
	}

	return commonVerts;
}

bool Parameterization::IsCutEdge(vector<int> cutVertices, vector<int> verts)
{
	if(verts.size() != 2)
		return false;
	
	for(int i = 0; i < cutVertices.size(); i++)
	{
		if(verts[0] == cutVertices[i])
		{
			if((i > 0 && verts[1] == cutVertices[i - 1]) || (i < cutVertices.size() - 1 && verts[1] == cutVertices[i + 1]))
				return true;
		}
	}

	return false;
}

pair<map<int, int>, set<int>> Parameterization::CreateCut(Mesh *mesh, vector<int> &cutVertices)
{
	// Store global Ids for cut vertices and triangles.
	FindCutVertices(mesh, cutVertices);
	vector<int> cutVertsWithoutStartEnd = cutVertices;
	cutVertsWithoutStartEnd.erase(cutVertsWithoutStartEnd.begin());
	cutVertsWithoutStartEnd.erase(cutVertsWithoutStartEnd.begin() + cutVertsWithoutStartEnd.size() - 1);
	set<int> cutTris;

	// Find all the triangles adjacent to the cut vertices except from the start and the end.
	for(int i = 1; i < cutVertices.size() - 1; i++)
	{
		vector<int> neighborTris = mesh->verts[cutVertices[i]]->triList;

		for(int j = 0; j < neighborTris.size(); j++)
		{
			cutTris.insert(neighborTris[j]);
		}
	}

	int label = 1;
	// Store (triId, triLabel) pairs
	map<int, int> triLabels;
	// Store triIds
	vector<int> consideredTris;
	vector<int> notConsideredTris(cutTris.begin(), cutTris.end());
		
	// Label the first triangle (cutTris[0])
	triLabels[*next(cutTris.begin(), 0)] = label;
	consideredTris.push_back(*next(cutTris.begin(), 0));

	// Consider all triangles around the cut.
	while(triLabels.size() != cutTris.size())
	{
		int cutTri = notConsideredTris[0];
		int triLabel;
		bool labelFound = false;

		for(int k = 0; k < consideredTris.size(); k++)
		{
			vector<int> commonVerts = FindCommonVertices(mesh, cutTri, consideredTris[k]);
			bool isCutEdge = IsCutEdge(cutVertices, commonVerts);
			if(isCutEdge == true && commonVerts.size() == 2)
			{
				// If common edge is a cut edge
				// Label different than consideredTris[k]
				triLabel = !triLabels[consideredTris[k]];
				labelFound = true;
				break;
			}
			else if(isCutEdge == false && commonVerts.size() == 2)
			{
				// If common edge is not a cut edge
				// Label same with consideredTris[k]
				triLabel = triLabels[consideredTris[k]];
				labelFound = true;
			}
		}

		if(labelFound == true)
		{
			triLabels[cutTri] = triLabel;
			consideredTris.push_back(cutTri);
			notConsideredTris.erase(notConsideredTris.begin());
		}
		else
		{
			notConsideredTris.erase(notConsideredTris.begin());
			notConsideredTris.push_back(cutTri);
		}
	}
	
	/***************If all the cut triangles are found and labeled, continue***************/
	// Store triangle Ids wrt labels.
	vector<int> labeled0;
	vector<int> labeled1;	
	FillLabels(triLabels, labeled0, labeled1);
	
	// Add duplicate cut vertices to mesh.
	vector<int> duplicateVertices;
	for(int i = 0; i < cutVertsWithoutStartEnd.size(); i++)
	{
		float *coords = mesh->verts[cutVertsWithoutStartEnd[i]]->coords;
		mesh->addVertex(coords[0], coords[1], coords[2]);
		int id = mesh->verts.size() - 1;
		duplicateVertices.push_back(id);
		// Update vertList.
		if(i == 0)		// First cut vertex
		{
			mesh->verts[id]->vertList.push_back(cutVertices[0]);
			mesh->verts[cutVertices[0]]->vertList.push_back(id);
		}
		else if(i == cutVertsWithoutStartEnd.size() - 1)	// Last cut vertex
		{
			mesh->verts[id]->vertList.push_back(id - 1);
			mesh->verts[id - 1]->vertList.push_back(id);
			
			mesh->verts[id]->vertList.push_back(cutVertices.back());
			mesh->verts[cutVertices.back()]->vertList.push_back(id);
		}
		else
		{
			mesh->verts[id]->vertList.push_back(id - 1);
			mesh->verts[id - 1]->vertList.push_back(id);
		}
	}

	cout << "number of duplicate vertices: " << duplicateVertices.size() << endl;
	
	// Arrange triangle vertex relationships.
	for(int i = 0; i < labeled0.size(); i++)
	{
		int triId = labeled0[i];
		int v1 = mesh->tris[labeled0[i]]->v1i;
		int v2 = mesh->tris[labeled0[i]]->v2i;
		int v3 = mesh->tris[labeled0[i]]->v3i;
		
		// Store cut and non-cut vertices.
		vector<int> cuts;	// without start and end
		vector<int> noncuts;
		FillCutsNonCuts(mesh, cuts, noncuts, cutVertices, triId);
		
	 	// duplicateVertices and cutVertsWithoutStartEnd have the same indexing.
		// Find index of v1
		int idx = GetIndex(cutVertsWithoutStartEnd, v1);
		if(idx > -1)
		{
			// Change cut vertex of 0-triangle from v1 to duplicate.
			mesh->tris[triId]->v1i = duplicateVertices[idx];
			
			for(int j = 0; j < noncuts.size(); j++)
			{
				// Add noncut of 0-triangle to duplicate.
				int noncutIdx = GetIndex(mesh->verts[duplicateVertices[idx]]->vertList, noncuts[j]);
				if(noncutIdx == -1)
					mesh->verts[duplicateVertices[idx]]->vertList.push_back(noncuts[j]);
				// Delete v1 from noncut of 0-triangle and add duplicate instead.
				int v1Idx = GetIndex(mesh->verts[noncuts[j]]->vertList, v1);
				if(v1Idx > -1)
				{
					mesh->verts[noncuts[j]]->vertList.erase(mesh->verts[noncuts[j]]->vertList.begin() + v1Idx);
					mesh->verts[noncuts[j]]->vertList.push_back(duplicateVertices[idx]);
				}
				// Delete noncut of 0-triangle from v1.
				noncutIdx = GetIndex(mesh->verts[v1]->vertList, noncuts[j]);
				if(noncutIdx > -1)
					mesh->verts[v1]->vertList.erase(mesh->verts[v1]->vertList.begin() + noncutIdx);
			}
			
			// Delete 0-triangle from v1 and add it to duplicate.
			int triIdx = GetIndex(mesh->verts[v1]->triList, triId);
			if(triIdx > -1)
			{
				mesh->verts[v1]->triList.erase(mesh->verts[v1]->triList.begin() + triIdx);
				mesh->verts[duplicateVertices[idx]]->triList.push_back(triId);
			}
		}
		// Find index of v2
		idx = GetIndex(cutVertsWithoutStartEnd, v2);
		if(idx > -1)
		{
			// Change cut vertex of 0-triangle from v2 to duplicate.
			mesh->tris[triId]->v2i = duplicateVertices[idx];
			
			for(int j = 0; j < noncuts.size(); j++)
			{
				// Add noncut of 0-triangle to duplicate.
				int noncutIdx = GetIndex(mesh->verts[duplicateVertices[idx]]->vertList, noncuts[j]);
				if(noncutIdx == -1)
					mesh->verts[duplicateVertices[idx]]->vertList.push_back(noncuts[j]);
				// Delete v2 from noncut of 0-triangle and add duplicate instead.
				int v2Idx = GetIndex(mesh->verts[noncuts[j]]->vertList, v2);
				if(v2Idx > -1)
				{
					mesh->verts[noncuts[j]]->vertList.erase(mesh->verts[noncuts[j]]->vertList.begin() + v2Idx);
					mesh->verts[noncuts[j]]->vertList.push_back(duplicateVertices[idx]);
				}
				// Delete noncut of 0-triangle from v2.
				noncutIdx = GetIndex(mesh->verts[v2]->vertList, noncuts[j]);
				if(noncutIdx > -1)
					mesh->verts[v2]->vertList.erase(mesh->verts[v2]->vertList.begin() + noncutIdx);
			}
			
			// Delete 0-triangle from v2 and add it to duplicate.
			int triIdx = GetIndex(mesh->verts[v2]->triList, triId);
			if(triIdx > -1)
			{
				mesh->verts[v2]->triList.erase(mesh->verts[v2]->triList.begin() + triIdx);
				mesh->verts[duplicateVertices[idx]]->triList.push_back(triId);
			}
		}
		// Find index of v3
		idx = GetIndex(cutVertsWithoutStartEnd, v3);
		if(idx > -1)
		{
			// Change cut vertex of 0-triangle from v3 to duplicate.
			mesh->tris[triId]->v3i = duplicateVertices[idx];
			
			for(int j = 0; j < noncuts.size(); j++)
			{
				// Add noncut of 0-triangle to duplicate.
				int noncutIdx = GetIndex(mesh->verts[duplicateVertices[idx]]->vertList, noncuts[j]);
				if(noncutIdx == -1)
					mesh->verts[duplicateVertices[idx]]->vertList.push_back(noncuts[j]);
				// Delete v3 from noncut of 0-triangle and add duplicate instead.
				int v3Idx = GetIndex(mesh->verts[noncuts[j]]->vertList, v3);
				if(v3Idx > -1)
				{
					mesh->verts[noncuts[j]]->vertList.erase(mesh->verts[noncuts[j]]->vertList.begin() + v3Idx);
					mesh->verts[noncuts[j]]->vertList.push_back(duplicateVertices[idx]);
				}
				// Delete noncut of 0-triangle from v3.
				noncutIdx = GetIndex(mesh->verts[v3]->vertList, noncuts[j]);
				if(noncutIdx > -1)
					mesh->verts[v3]->vertList.erase(mesh->verts[v3]->vertList.begin() + noncutIdx);
			}
			
			// Delete 0-triangle from v3 and add it to duplicate.
			int triIdx = GetIndex(mesh->verts[v3]->triList, triId);
			if(triIdx > -1)
			{
				mesh->verts[v3]->triList.erase(mesh->verts[v3]->triList.begin() + triIdx);
				mesh->verts[duplicateVertices[idx]]->triList.push_back(triId);
			}
		}
	}
	
	// Add duplicateVertices to cutVertices to find the whole closed boundary.
	reverse(duplicateVertices.begin(), duplicateVertices.end());
	cutVertices.insert(cutVertices.end(), duplicateVertices.begin(), duplicateVertices.end());

	return make_pair(triLabels, cutTris);
}

void Parameterization::FillCutsNonCuts(Mesh *mesh, vector<int> &cuts, vector<int> &noncuts, vector<int> cutVertices, int triId)
{
	vector<int> triVerts;
	triVerts.push_back(mesh->tris[triId]->v1i);
	triVerts.push_back(mesh->tris[triId]->v2i);
	triVerts.push_back(mesh->tris[triId]->v3i);
	
	for(int i = 0; i < triVerts.size(); i++)
	{
		int idx = GetIndex(cutVertices, triVerts[i]);
		if(idx == 0 || idx == cutVertices.size() - 1)	// If start or end vertex
		{
			// Skip
		}
		else if(idx > -1)
		{
			cuts.push_back(triVerts[i]);
		}
		else
		{
			noncuts.push_back(triVerts[i]);
		}
	}
}
		    
void Parameterization::FillLabels(map<int, int> triLabels, vector<int> &labeled0, vector<int> &labeled1)
{
	map<int, int>::iterator it = triLabels.begin();
	for(; it != triLabels.end(); it++)
	{
		if(it->second == 0)
		{
			labeled0.push_back(it->first);
		}
		else	
		{
			labeled1.push_back(it->first);
		}
	}
}

vector<FreeVertex *> Parameterization::DiskParameterization(Mesh *mesh)
{
	vector<FreeVertex *> resultVertices;

	// Get boundary vertices
	set<int> boundaryVertices = FindBoundaryVertices(mesh);
	vector<int> boundaryVerts(boundaryVertices.begin(), boundaryVertices.end());
	// Map boundary vertices to disk
	vector<FreeVertex *> circleVertices = GetCircleVertices(boundaryVertices.size());

	// Map non-boundary vertices
	if(weight == UNIFORM)
	{
		resultVertices = ParamUniform(mesh, circleVertices, boundaryVerts);
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

vector<FreeVertex *> Parameterization::ParamClosedMesh(Mesh *mesh, vector<int> boundaryVertices)
{
	vector<FreeVertex *> resultVertices;
	// Order does not matter for SetBoundaryVertices function.
	set<int> boundaryVerts(boundaryVertices.begin(), boundaryVertices.end());
	SetBoundaryVertices(mesh, boundaryVerts);

	// Map boundary vertices to disk
	vector<FreeVertex *> circleVertices = GetCircleVertices(boundaryVertices.size());

	// Map non-boundary vertices
	if(weight == UNIFORM)
	{
		resultVertices = ParamUniform(mesh, circleVertices, boundaryVertices);
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
