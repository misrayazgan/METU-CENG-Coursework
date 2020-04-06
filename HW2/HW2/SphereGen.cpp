#include "SphereGen.h"

FreeVertex* SphereGen::FindCentroid(Mesh *mesh)
{
	float *coords = new float[3];
	coords[0] = 0;
	coords[1] = 0;
	coords[2] = 0;

	for(int i = 0; i < mesh->verts.size(); i++)
	{
		coords[0] += mesh->verts[i]->coords[0];
		coords[1] += mesh->verts[i]->coords[1];
		coords[2] += mesh->verts[i]->coords[2];
	}

	coords[0] /= mesh->verts.size();
	coords[1] /= mesh->verts.size();
	coords[2] /= mesh->verts.size();

	return new FreeVertex(coords);
}

vector<int> SphereGen::FindPoleVertices(Mesh *mesh)
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

void SphereGen::FindCutVertices(Mesh *mesh, vector<int> &cutVertices)
{
	// Find start and end points of the cut.
	vector<int> poles = FindPoleVertices(mesh);
	Dijkstra *d = new Dijkstra();
	cutVertices = d->GetShortestPath(mesh, poles[0], poles[1]);
}

// Given two triangles find the common edge vertices.
vector<int> SphereGen::FindCommonVertices(Mesh *mesh, int t1, int t2)
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

bool SphereGen::IsCutEdge(vector<int> cutVertices, vector<int> verts)
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

// If v is in cutVertices, return its index.
// Else return -1.
int SphereGen::GetIndex(const vector<int> cutVertices, int v)
{
	//ptrdiff_t
	int result = find(cutVertices.begin(), cutVertices.end(), v) - cutVertices.begin();
	if(result < cutVertices.size())
	{
		return result; 
	}
	return -1;
}

pair<map<int, int>, set<int>> SphereGen::CreateCut(Mesh *mesh, vector<int> &cutVertices, vector<int> &duplicateVertices)
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
	//for(int i = 1; i < cutTris.size(); i++)
	{
		int cutTri = notConsideredTris[0];//*next(cutTris.begin(), i);
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
	
	// Add duplicate cut vertices.
	//vector<int> duplicateVertices;
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

void SphereGen::FillCutsNonCuts(Mesh *mesh, vector<int> &cuts, vector<int> &noncuts, vector<int> cutVertices, int triId)
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
		    
void SphereGen::FillLabels(map<int, int> triLabels, vector<int> &labeled0, vector<int> &labeled1)
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

float * SphereGen::GetMiddlePoint(float *v1coords, float *v2coords)
{
	float *midCoords = new float[3];
	midCoords[0] = (v1coords[0] + v2coords[0]) / 2;
	midCoords[1] = (v1coords[1] + v2coords[1]) / 2;
	midCoords[2] = (v1coords[2] + v2coords[2]) / 2;

	return midCoords;
}


void SphereGen::DivideTriangles(Mesh *mesh, int v1, int v2, int v3, int level)
{
	if(level-- > 0)
	{
		// Find coordinates for middle of 1 and 2
		float *mid12 = GetMiddlePoint(mesh->verts[v1]->coords, mesh->verts[v2]->coords);
		int id12 = mesh->findVertexByCoord(mid12);
		if(id12 > -1)	// vertex is already added to mesh
		{
			// Continue with id12
		}
		else		// vertex is not added to mesh previously
		{
			mesh->addVertex(mid12[0], mid12[1], mid12[2]);
			id12 = mesh->verts.size() - 1;
		}

		// Find coordinates for middle of 2 and 3
		float *mid23 = GetMiddlePoint(mesh->verts[v2]->coords, mesh->verts[v3]->coords);
		int id23 = mesh->findVertexByCoord(mid23);
		if(id23 > -1)	// vertex is already added to mesh
		{
			// Continue with id23
		}
		else		// vertex is not added to mesh previously
		{
			mesh->addVertex(mid23[0], mid23[1], mid23[2]);
			id23 = mesh->verts.size() - 1;
		}

		// Find coordinates for middle of 1 and 3
		float *mid13 = GetMiddlePoint(mesh->verts[v1]->coords, mesh->verts[v3]->coords);
		int id13 = mesh->findVertexByCoord(mid13);
		if(id13 > -1)	// vertex is already added to mesh
		{
			// Continue with id13
		}
		else		// vertex is not added to mesh previously
		{
			mesh->addVertex(mid13[0], mid13[1], mid13[2]);
			id13 = mesh->verts.size() - 1;
		}

		DivideTriangles(mesh, v1, id12, id13, level);
		DivideTriangles(mesh, id12, v2, id23, level);
		DivideTriangles(mesh, id12, id23, id13, level);
		DivideTriangles(mesh, id13, id23, v3, level);
	}
	else
	{
		mesh->addTriangle(v1, v2, v3);
	}
}

void SphereGen::SplitTetrahedron(Mesh *mesh, int level)
{
	DivideTriangles(mesh, 0, 1, 2, level);
	DivideTriangles(mesh, 0, 2, 3, level);
	DivideTriangles(mesh, 0, 3, 1, level);
	DivideTriangles(mesh, 1, 3, 2, level);
}

Mesh * SphereGen::CreateTetrahedron()
{
	Mesh *mesh = new Mesh();
	float c0[3] = {0, 0, 0};
	float c1[3] = {1, 0, 0};
	float c2[3] = {0, 1, 0};
	float c3[3] = {0, 0, 1};
	float *coords0 = c0;
	float *coords1 = c1;
	float *coords2 = c2;
	float *coords3 = c3;
	
	mesh->addVertex(coords0[0], coords0[1], coords0[2]);
	mesh->addVertex(coords1[0], coords1[1], coords1[2]);
	mesh->addVertex(coords2[0], coords2[1], coords2[2]);
	mesh->addVertex(coords3[0], coords3[1], coords3[2]);

	return mesh;
}

void SphereGen::GetSphereVertices(Mesh *mesh)
{
	float radius = 1;

	for(int i = 0; i < mesh->verts.size(); i++)
	{
		Vertex *v = mesh->verts[i];
		// Compute the vector
		float *sphereCoords = new float[3];
		sphereCoords[0] = v->coords[0] - centroid->coords[0];
		sphereCoords[1] = v->coords[1] - centroid->coords[1];
		sphereCoords[2] = v->coords[2] - centroid->coords[2];

		// Normalize the vector
		float len = sqrt(pow(sphereCoords[0], 2) + pow(sphereCoords[1], 2) + pow(sphereCoords[2], 2));
		sphereCoords[0] /= len;
		sphereCoords[1] /= len;
		sphereCoords[2] /= len;

		// Scale wrt radius
		sphereCoords[0] *= radius;
		sphereCoords[1] *= radius;
		sphereCoords[2] *= radius;

		mesh->verts[i]->coords = sphereCoords;
	}
}

Mesh * SphereGen::GenerateSphere()
{
	// Start with 4 vertices of the tetrahedron.
	Mesh *mesh = CreateTetrahedron();
	// Split each face to 4 triangles recursively.
	int level = 5;
	SplitTetrahedron(mesh, level);
	// Find the center point of the tetrahedron.
	centroid = FindCentroid(mesh);

	GetSphereVertices(mesh);
	return mesh;
}
