#include "ConvexHull.h"
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"

void ConvexHull::FindConvexHull(Mesh *mesh)
{
	int numberOfVertices = mesh->verts.size();
	ch_vertex* vertices = new ch_vertex[numberOfVertices];

	for(int i = 0; i < mesh->verts.size(); i++)
	{
		vertices[i].x = mesh->verts[i]->coords[0];
		vertices[i].y = mesh->verts[i]->coords[1];
		vertices[i].z = mesh->verts[i]->coords[2];
	}

	int* faceIndices = NULL;
	int nFaces;
	convhull_3d_build(vertices, numberOfVertices, &faceIndices, &nFaces);
	/* Where 'faceIndices' is a flat 2D matrix [nFaces x 3] */
	char *name = "test";
	convhull_3d_export_obj(vertices, numberOfVertices, faceIndices, nFaces, 1, name);
	free(vertices);
	free(faceIndices);
}

//#define QUICKHULL_IMPLEMENTATION
//#include "quickhull.h"
//
//
//void ConvexHull::FindConvexHull(Mesh *mesh)
//{
//	const int a = mesh->verts.size();
//	qh_vertex_t *vertices = new qh_vertex_t[a];
//
//	for(int i = 0; i < mesh->verts.size(); i++)
//	{
//		vertices[i].x = mesh->verts[i]->coords[0];
//		vertices[i].y = mesh->verts[i]->coords[1];
//		vertices[i].z = mesh->verts[i]->coords[2];
//	}
//
//	qh_mesh_t hull_mesh = qh_quickhull3d(vertices, mesh->verts.size());
//
//	qh_mesh_export(&hull_mesh, "test.obj");
//	qh_free_mesh(hull_mesh);
//}

