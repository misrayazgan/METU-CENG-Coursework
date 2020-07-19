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