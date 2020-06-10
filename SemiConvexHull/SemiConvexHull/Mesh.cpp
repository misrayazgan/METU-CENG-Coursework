#include "Mesh.h"


void Mesh::loadOff(char* name)
{
	FILE* fPtr = fopen(name, "r");
	char str[334];

	fscanf(fPtr, "%s", str);

	int nVerts, nTris, n, i = 0;
	float x, y, z;

	fscanf(fPtr, "%d %d %d\n", &nVerts, &nTris, &n);
	while (i++ < nVerts)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addVertex(x, y, z);
	}

	while (fscanf(fPtr, "%d", &i) != EOF)
	{
		fscanf(fPtr, "%f %f %f", &x, &y, &z);
		addTriangle((int) x, (int) y, (int) z);
	}

	computeNormals();

	fclose(fPtr);
}


void Mesh::loadObj(char* name)
{
	ifstream file (name);
	string line;
	string lineId;
	float x,y,z;
	string xs,ys,zs;
	string xt, yt, zt;
	vector<int> vertIds;	// store vertexIds line by line

	while (!file.eof())
	{
		getline(file, line);
		
		if (line[0] == 'v' && line[1] == ' ')
		{
			// Vertices
			std::istringstream iss(line);
			iss >> lineId >> x >> y >> z;

			// Check if the vertex is already added.
			float *coords = new float[3];
			coords[0] = x;
			coords[1] = y;
			coords[2] = z;
			int v = findVertexByCoord(coords);
			if(v > -1)
			{
				// Vertex was already added.
				vertIds.push_back(v);
			}
			else
			{
				// Vertex was not found, add new vertex.
				addVertex(x, y, z);
				vertIds.push_back(verts.size() - 1);
			}
		}
		else if(line[0] == 'f')
		{
			// Triangles
			std::istringstream iss(line);
			iss >> lineId >> xt >> yt >> zt;
			parseString(xt, xs);
			parseString(yt, ys);
			parseString(zt, zs);
			addTriangle(vertIds[stoi(xs) - 1], vertIds[stoi(ys) - 1], vertIds[stoi(zs) - 1]);
		}
	}

	computeNormals();
}


void Mesh::parseString(string &src, string &dest)
{
	for(int i = 0; i < src.size(); i++)
	{
		if(src[i] == '/')
		{
			dest.assign(src, 0, i);
			break;
		}
	}
}

void Mesh::createCube(float sideLen)
{
	//coordinates
	float flbc[3] = {0, 0, 0}, deltaX = 0, deltaY = 0, deltaZ = 0;
	for (int v = 0; v < 8; v++)
	{
		switch (v)
		{
			case 1:
				deltaX = sideLen;
				break;
			case 2:
				deltaZ = -sideLen;
				break;
			case 3:
				deltaX = 0;
				break;
			case 4:
				deltaZ = 0;
				deltaY = sideLen;
				break;
			case 5:
				deltaX = sideLen;
				break;
			case 6:
				deltaZ = -sideLen;
				break;
			default:
				deltaX = 0;;
				break;
		}
		addVertex(flbc[0] + deltaX, flbc[1] + deltaY, flbc[2] + deltaZ);
	}

	addTriangle(0, 2, 1);
	addTriangle(0, 3, 2);

	addTriangle(1, 2, 5);
	addTriangle(2, 6, 5);

	addTriangle(2, 3, 6);
	addTriangle(3, 7, 6);

	addTriangle(3, 4, 7);
	addTriangle(3, 0, 4);

	addTriangle(4, 5, 6);
	addTriangle(4, 6, 7);

	addTriangle(0, 1, 5);
	addTriangle(0, 5, 4);
}

void Mesh::addTriangle(int v1, int v2, int v3)
{
	int idx = tris.size();
	tris.push_back( new Triangle(idx, v1, v2, v3) );

	//set up structure

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);

	if (! makeVertsNeighbor(v1, v2) )
		addEdge(v1, v2);

	if (! makeVertsNeighbor(v1, v3) )
		addEdge(v1, v3);

	if (! makeVertsNeighbor(v2, v3) )
		addEdge(v2, v3);

}

void Mesh::addTriangleWithoutEdge(int v1, int v2, int v3)
{
	int idx = tris.size();
	tris.push_back( new Triangle(idx, v1, v2, v3) );

	verts[v1]->triList.push_back(idx);
	verts[v2]->triList.push_back(idx);
	verts[v3]->triList.push_back(idx);
}

bool Mesh::makeVertsNeighbor(int v1i, int v2i)
{
	//returns true if v1i already neighbor w/ v2i; false o/w

	for (int i = 0; i < verts[v1i]->vertList.size(); i++)
		if (verts[v1i]->vertList[i] == v2i)
			return true;


	verts[v1i]->vertList.push_back(v2i);
	verts[v2i]->vertList.push_back(v1i);
	return false;
}

void Mesh::addVertex(float x, float y, float z)
{
	int idx = verts.size();
	float* c = new float[3];
	c[0] = x;
	c[1] = y;
	c[2] = z;

	verts.push_back( new Vertex(idx, c) );
}

void Mesh::addEdge(int v1, int v2)
{
	int idx = edges.size();

	float len = computeEdgeLength(v1, v2);
	edges.push_back(new Edge(idx, v1, v2, len));

	verts[v1]->edgeList.push_back(idx);
	verts[v2]->edgeList.push_back(idx);
}

float Mesh::computeEdgeLength(int v1i, int v2i)
{
	Vertex *v1 = verts[v1i];
	Vertex *v2 = verts[v2i];

	float *v1Coords = v1->coords;
	float *v2Coords = v2->coords;

	float len = sqrt(pow((v2Coords[0] - v1Coords[0]), 2) + pow((v2Coords[1] - v1Coords[1]), 2) + pow((v2Coords[2] - v1Coords[2]), 2));
	return len;
}

bool Mesh::isNeighbor(int v1i, int v2i)
{
	for (int i = 0; i < verts[v1i]->vertList.size(); i++)
		if (verts[v1i]->vertList[i] == v2i)
			return true;

	return false;
}

float * Mesh::crossProduct(float *a, float *b)
{
    float *result = new float[3];
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];

    return result;
}

void Mesh::computeNormals()
{
	for(int i = 0; i < tris.size(); i++)
	{
		float *v1Coords = verts[tris[i]->v1i]->coords;
		float *v2Coords = verts[tris[i]->v2i]->coords;
		float *v3Coords = verts[tris[i]->v3i]->coords;

		// b - a
		float *bminusa = new float[3];
		bminusa[0] = v2Coords[0] - v1Coords[0];
		bminusa[1] = v2Coords[1] - v1Coords[1];
		bminusa[2] = v2Coords[2] - v1Coords[2];

		// c - a
		float *cminusa = new float[3];
		cminusa[0] = v3Coords[0] - v1Coords[0];
		cminusa[1] = v3Coords[1] - v1Coords[1];
		cminusa[2] = v3Coords[2] - v1Coords[2];

		// cross product
		float *normal = crossProduct(bminusa, cminusa);

		// normalize
		float len = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
		normal[0] /= len;
		normal[1] /= len;
		normal[2] /= len;

		tris[i]->normal = normal;
	}
}

int Mesh::findVertexByCoord(float *coords)
{
	for(int i = 0; i < verts.size(); i++)
	{
		float *vCoords = verts[i]->coords;
		if(coords[0] == vCoords[0] && coords[1] == vCoords[1] && coords[2] == vCoords[2])
		{
			return i;
		}
	}
	return -1;
}

int Mesh::findEdgeByVertices(int v1, int v2)
{
	for(int i = 0; i < edges.size(); i++)
	{
		if((edges[i]->v1i == v1 && edges[i]->v2i == v2) ||
			(edges[i]->v1i == v2 && edges[i]->v2i == v1))
		{
			return i;
		}
	}

	return -1;
}

vector<int> Mesh::findCommonTriangles(int v1, int v2)
{
	vector<int> commonTris;
	for(int i = 0; i < verts[v1]->triList.size(); i++)
	{
		int tri = verts[v1]->triList[i];

		for(int j = 0; j < verts[v2]->triList.size(); j++)
		{
			if(tri == verts[v2]->triList[j])
			{
				commonTris.push_back(tri);
				break;
			}
		}
	}

	return commonTris;
}

int Mesh::getThirdVertexOfTriangle(int t, int v1, int v2)
{
	Triangle *tri = tris[t];
	if(tri->v1i != v1 && tri->v1i != v2)
		return tri->v1i;
	if(tri->v2i != v1 && tri->v2i != v2)
		return tri->v2i;
	if(tri->v3i != v1 && tri->v3i != v2)
		return tri->v3i;
}

int Mesh::getTriangleByVertices(int v1, int v2, int v3)
{
	vector<int> vertices;
	vertices.push_back(v1);
	vertices.push_back(v2);
	vertices.push_back(v3);

	sort(vertices.begin(), vertices.end());

	for(int i = 0; i < tris.size(); i++)
	{
		vector<int> triVerts;
		triVerts.push_back(tris[i]->v1i);
		triVerts.push_back(tris[i]->v2i);
		triVerts.push_back(tris[i]->v3i);

		sort(triVerts.begin(), triVerts.end());
		
		if(vertices == triVerts)
		{
			return i;
		}
	}

	return -1;
}

float Mesh::computeAvgEdgeLength()
{
	float sum = 0;
	for(int i = 0; i < edges.size(); i++)
	{
		sum += edges[i]->length;
	}

	sum /= edges.size();
	return sum;
}
