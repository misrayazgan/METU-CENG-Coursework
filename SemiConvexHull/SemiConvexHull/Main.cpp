#include "Remeshing.h"
#include "ConvexHull.h"
#include "Sampling.h"
#include "Painter.h"

#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

int main(int, char ** argv)
{
	HWND window = SoWin::init(argv[0]);
	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);
	SoSeparator * root = new SoSeparator;
	root->ref();

	float * blue = new float[3];
	blue[0] = 0.0f;
	blue[1] = 0.0f;
	blue[2] = 1.0f;

	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();

	char * filename = "horse0.off";
	mesh->loadOff(filename);

	cout << "Number of Vertices: " << mesh->verts.size() << endl;
	cout << "Number of Edges: " << mesh->edges.size() << endl;
	cout << "Number of Triangles: " << mesh->tris.size() << endl;

	float avgLen = mesh->computeAvgEdgeLength();

	/******************Sample points on the mesh surface uniformly******************/
	/*root->addChild(painter->getShapeSep(mesh));

	Sampling *sampling = new Sampling(mesh->tris.size());
	vector<float *> samples = sampling->UniformSampling(mesh);

	for(int i = 0; i < samples.size(); i++)
	{
		root->addChild(painter->getSphereSepByCoord(mesh, samples[i], 0.5f));
	}*/

	/******************Find convex hull of the mesh*******************/
	ConvexHull *hull = new ConvexHull();
	hull->FindConvexHull(mesh);

	Mesh *hullMesh = new Mesh();
	hullMesh->loadObj("test.obj");
	//hullMesh = mesh;
	//hullMesh->loadOff("square.off");
	//hullMesh->createCube(1.0f);
	
	cout << "Hull mesh Number of Vertices: " << hullMesh->verts.size() << endl;
	cout << "Hull mesh Number of Edges: " << hullMesh->edges.size() << endl;
	cout << "Hull mesh Number of Triangles: " << hullMesh->tris.size() << endl;

	//root->addChild(painter->getShapeSep(hullMesh));

	/*************************CONTROL FOR CONVEX HULL******************************/
	vector<vector<int>> edgeCount(hullMesh->verts.size(), vector<int>(hullMesh->verts.size(), 0));

	for(int k = 0; k < hullMesh->tris.size(); k++)
	{
		Triangle *tri = hullMesh->tris[k];
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

	for(int i = 0; i < edgeCount.size(); i++)
	{
		for(int j = 0; j < edgeCount.size(); j++)
		{
			if(!(edgeCount[i][j] == 0 || edgeCount[i][j] == 2))
			{
				int a = 0;
			}
		}
	}
	/*************************CONTROL FOR CONVEX HULL******************************/
	/*root->addChild(painter->getShapeSep(hullMesh));
	for(int i = 0; i < hullMesh->verts.size(); i++)
	{
		root->addChild(painter->getSphereSep(hullMesh, i, 1, 0.2f));
	}
	for(int i = 0; i < hullMesh->edges.size(); i++)
	{
		float * v1coords = hullMesh->verts[hullMesh->edges[i]->v1i]->coords;
		float * v2coords = hullMesh->verts[hullMesh->edges[i]->v2i]->coords;
		root->addChild(painter->drawEdge(hullMesh, v1coords, v2coords, blue, false));
	}*/

	/*****************Remeshing the convex hull******************/
	Sampling *sampling = new Sampling(hullMesh->tris.size() * 5);
	vector<float *> samples = sampling->UniformSampling(hullMesh);

	Remeshing *remesher = new Remeshing();
	remesher->Remesh(hullMesh, avgLen);
	remesher->OptimizeMesh(hullMesh, avgLen, samples);

	cout << "AFTER REMESHING" << endl;
	cout << "Hull mesh Number of Vertices: " << hullMesh->verts.size() << endl;
	cout << "Hull mesh Number of Edges: " << hullMesh->edges.size() << endl;
	cout << "Hull mesh Number of Triangles: " << hullMesh->tris.size() << endl;

	root->addChild(painter->getShapeSep(hullMesh));
	for(int i = 0; i < hullMesh->edges.size(); i++)
	{
		float * v1coords = hullMesh->verts[hullMesh->edges[i]->v1i]->coords;
		float * v2coords = hullMesh->verts[hullMesh->edges[i]->v2i]->coords;
		root->addChild(painter->drawEdge(hullMesh, v1coords, v2coords, blue, true));
	}






	/*****************Gradient Descent***************/



	viewer->setSize(SbVec2s(1000, 1000));
	//viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}