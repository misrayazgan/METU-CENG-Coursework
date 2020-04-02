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

vector<int> SphereGen::FindCutVertices(Mesh *mesh)
{
	// Find start and end points of the cut.
	vector<int> poles = FindPoleVertices(mesh);
	Dijkstra *d = new Dijkstra();
	vector<int> cutVertices = d->GetShortestPath(mesh, poles[0], poles[1]);
	return cutVertices;
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

pair<vector<int>, set<int>> SphereGen::CreateCut(Mesh *mesh)
{
	// Store global Ids for cut vertices and triangles.
	vector<int> cutVertices = FindCutVertices(mesh);
	vector<int> cutVertsWithoutStartEnd = cutVertices;
	cutVertsWithoutStartEnd.erase(cutVertsWithoutStartEnd.begin());
	cutVertsWithoutStartEnd.erase(cutVertsWithoutStartEnd.end());
	set<int> cutTris;

	// Duplicate cut vertices except from the start and the end.
	// Find all the triangles adjacent to the cut vertices.
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
	vector<pair<int, int>> triLabels;
	vector<int> consideredTris;
		
	// Label the first triangle (cutTris[0])
	triLabels.push_back(make_pair(0, label));
	consideredTris.push_back(*next(cutTris.begin(), 0));

	// Consider all triangles around the cut.
	for(int i = 1; i < cutTris.size(); i++)
	{
		int cutTri = *next(cutTris.begin(), i);
		int triLabel;
		bool labelFound = false;

		for(int k = 0; k < consideredTris.size(); k++)
		{
			vector<int> commonVerts = FindCommonVertices(mesh, cutTri, consideredTris[k]);
			bool isCutEdge = IsCutEdge(cutVertices, commonVerts);
			if(isCutEdge == true)
			{
				// If common edge is a cut edge
				// Label different than consideredTris[k]
				triLabel = !triLabels[k];
				labelFound = true;
				break;
			}
			else if(isCutEdge == false && commonVerts.size() == 2)
			{
				// If common edge is not a cut edge
				// Label same with consideredTris[k]
				triLabel = triLabels[k];
				labelFound = true;
			}
			else
			{
				// If triangles have no edge in common(maybe common vertex).
				// Cannot decide the label, skip the triangle.
				labelFound = false;
			}
		}

		if(labelFound = true)
		{
			triLabels.push_back(make_pair(cutTri, triLabel));
		}
	}
	
	/***************to be deleted**************/
	vector<int> notLabeled;
	// Check if all cutTriIds are labeled. (might just check sizes of trilabels and cuttris)
	for(int i = 0; i < triLabels.size(); i++)
	{
		// Wrong but works for now
		if(triLabels[i].first != cutTris[i])
		{
			notLabeled.push_back(cutTris[i]);
		}
	}
	
	cout << "number of cut triangles: " << cutTris.size() << endl;
	cout << "number of not labeled cut triangles: " << notLabeled.size() << endl;
	/***************to be deleted**************/
	
	/***************If all the cut triangles are found and labeled, continue***************/
	// Store triangle Ids wrt labels.
	vector<int> labeled0;
	vector<int> labeled1;	
	FillLabels(triLabels, labeled0, labeled1);
	
	// Add duplicate cut vertices.
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
	
	// Arrange triangle vertex relationships.
	for(int i = 0; i < labeled0.size(); i++)
	{
		int triId = labeled0[i];
		int v1 = mesh->tris[labeled0[i]]->v1i;
		int v2 = mesh->tris[labeled0[i]]->v2i;
		int v3 = mesh->tris[labeled0[i]]->v3i;
		
		// Store cut and non-cut vertices.
		vector<int> cuts;
		vector<int> noncuts;
		FillCutsNonCuts(cuts, noncuts, cutVertices, triId);
		
		// Find index of v1
		int idx = GetIndex(cutVertsWithoutStartEnd, v1);
		if(idx > -1)
		{
			// duplicateVertices and cutVertsWithoutStartEnd have the same indexing.
			// Update vertices of triangle.
			mesh->tris[triId]->v1i = idx;
			
			for(int j = 0; j < noncuts.size(); j++)
			{
				// Update neighboring noncut vertices of cut vertices.
				mesh->verts[duplicateVertices[idx]]->vertList.push_back(noncuts[j]);
				// Update neighboring vertices of non-cut vertices.
				mesh->verts[noncuts[j]]->vertList.erase(v1);
			}
			
			// Update neighboring triangles of cut vertices.
			// erase element ??????????????????????
			mesh->verts[v1]->triList.erase(i);			
			mesh->verts[duplicateVertices[idx]]->triList.push_back(i);			
		}
		// Find index of v2
		int idx = GetIndex(cutVertsWithoutStartEnd, v2);
		if(idx > -1)
		{
			// duplicateVertices and cutVertsWithoutStartEnd have the same indexing.
			// Update vertices of triangle.
			mesh->tris[triId]->v2i = idx;
			
			for(int j = 0; j < noncuts.size(); j++)
			{
				// Update neighboring noncut vertices of cut vertices.
				mesh->verts[duplicateVertices[idx]]->vertList.push_back(noncuts[j]);
				// Update neighboring vertices of non-cut vertices.
				mesh->verts[noncuts[j]]->vertList.erase(v2);
			}
			
			// Update neighboring triangles of cut vertices.
			// erase element ??????????????????????
			mesh->verts[v2]->triList.erase(i);			
			mesh->verts[duplicateVertices[idx]]->triList.push_back(i);			
		}
		// Find index of v3
		int idx = GetIndex(cutVertsWithoutStartEnd, v3);
		if(idx > -1)
		{
			// duplicateVertices and cutVertsWithoutStartEnd have the same indexing.
			// Update vertices of triangle.
			mesh->tris[triId]->v3i = idx;
			
			for(int j = 0; j < noncuts.size(); j++)
			{
				// Update neighboring noncut vertices of cut vertices.
				mesh->verts[duplicateVertices[idx]]->vertList.push_back(noncuts[j]);
				// Update neighboring vertices of non-cut vertices.
				mesh->verts[noncuts[j]]->vertList.erase(v3);
			}
			
			// Update neighboring triangles of cut vertices.
			// erase element ??????????????????????
			mesh->verts[v3]->triList.erase(i);			
			mesh->verts[duplicateVertices[idx]]->triList.push_back(i);			
		}
	}	

	return make_pair(triLabels, cutTris);
}

void SphereGen::FillCutsNonCuts(vector<int> &cuts, vector<int> &noncuts, vector<int> cutVertices, int triId)
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
		    
void SphereGen::FillLabels(pair<vector<int>, vector<int>> triLabels, vector<int> &labeled0, vector<int> &labeled1)
{
	for(int i = 0; i < triLabels.size(); i++)
	{
		if(triLabels[i].second == 0)
		{
			labeled0.push_back(triLabels[i].first);
		}
		else	
		{
			labeled1.push_back(triLabels[i].first);
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
	if(level == 0)
	{
		mesh->addTriangle(v1, v2, v3);
	}
	else
	{
		// Find coordinates for middle of 1 and 2
		float *mid12 = GetMiddlePoint(mesh->verts[v1]->coords, mesh->verts[v2]->coords);
		mesh->addVertex(mid12[0], mid12[1], mid12[2]);
		int id12 = mesh->verts.size() - 1;

		// Find coordinates for middle of 1 and 2
		float *mid23 = GetMiddlePoint(mesh->verts[v2]->coords, mesh->verts[v3]->coords);
		mesh->addVertex(mid23[0], mid23[1], mid23[2]);
		int id23 = mesh->verts.size() - 1;

		// Find coordinates for middle of 1 and 2
		float *mid13 = GetMiddlePoint(mesh->verts[v1]->coords, mesh->verts[v3]->coords);
		mesh->addVertex(mid13[0], mid13[1], mid13[2]);
		int id13 = mesh->verts.size() - 1;

		DivideTriangles(mesh, v1, id12, id13, level - 1);
		DivideTriangles(mesh, id12, v2, id23, level - 1);
		DivideTriangles(mesh, id12, id23, id13, level - 1);
		DivideTriangles(mesh, id13, id23, v3, level - 1);

		level--;
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

