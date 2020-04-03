#pragma once

#include <iostream>
#include <vector>

using namespace std;

struct Vertex
{
	float* coords, * normals; //3d coordinates etc
	int idx; //who am i; verts[idx]

	vector< int > vertList; //adj vvertices;
	vector< int > triList; 
	vector< int > edgeList; 
	
	Vertex(int i, float* c) : idx(i), coords(c) {};
};

struct Edge
{
	int idx; //edges[idx]
	int v1i, v2i; //endpnts
	float length;
	Edge(int id, int v1, int v2, float len) : idx(id), v1i(v1), v2i(v2), length(len) { /*computeLength();*/ };

	/* void computeLength()
	{
		length = 7;
	}*/
};

struct Triangle
{
	int idx; //tris[idx]
	int v1i, v2i, v3i;
	float *normal;	//normal vector of face
	Triangle(int id, int v1, int v2, int v3) : idx(id), v1i(v1), v2i(v2), v3i(v3) {};
};

class Mesh
{
private:
	bool makeVertsNeighbor(int v1i, int v2i);
	float * crossProduct(float *a, float *b);
public:
	vector< Vertex* > verts;
	vector< Triangle* > tris;
	vector< Edge* > edges;


	Mesh() {} ;
	void addTriangle(int v1, int v2, int v3);
	void addEdge(int v1, int v2);
	void addVertex(float x, float y, float z);
	void createCube(float side);
	void loadOff(char* name);
	float computeEdgeLength(int v1i, int v2i);
	bool isNeighbor(int v1i, int v2i);
	void computeNormals();
	int findVertexByCoord(float *coords);
};

// Vertex without an idx, independent from mesh
struct FreeVertex
{
	float* coords; // 3d coordinates
	FreeVertex(float* c) : coords(c) {};
};

// Edge without an idx, independent from mesh
struct FreeEdge
{
	int v1i, v2i; //endpnts
	float length;
	FreeEdge(int v1, int v2, float len) : v1i(v1), v2i(v2), length(len) {};
};

typedef struct Vec3f
{
	float x;
	float y;
	float z;
};
