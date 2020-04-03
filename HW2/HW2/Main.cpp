#include "Parameterization.h"
#include <iomanip>
#include "SphereGen.h"
#include "SphereParam.h"
#include "Painter.h"

#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>


int main(int, char ** argv)
{
	HWND window = SoWin::init(argv[0]);

	SoWinExaminerViewer *viewer = new SoWinExaminerViewer(window);

	//make a dead simple scene graph by using the Coin library, only containing a single cone under the scenegraph root
	SoSeparator *root = new SoSeparator;
	root->ref();

	int taskId;
	cout << "1: Disk Parameterization | Uniform Weights" << endl;
	cout << "2: Disk Parameterization | Harmonic Weights" << endl;
	cout << "3: Disk Parameterization | Mean Value Weights" << endl;
	cout << "4: Spherical Parameterization of Closed Mesh" << endl;
	cout << "5: Disk Parameterization of Closed Mesh" << endl;
	cout << "6: Sphere Generation" << endl;
	cout << "7: Disk Parameterization of Generated Sphere" << endl;
	cout << "Enter the id for the desired task:" << endl;
	cin >> taskId;

	Mesh *mesh = new Mesh();
	char filename[20];
	if(taskId != 6 && taskId != 7)
	{
		cout << "Enter filename:" << endl;
		cin >> setw(20) >> filename;
		mesh->loadOff(filename);
	}

	

	cout << "Number of Vertices: " << mesh->verts.size() << endl;
	cout << "Number of Edges: " << mesh->edges.size() << endl;
	cout << "Number of Triangles: " << mesh->tris.size() << endl;

	Painter *painter = new Painter();
	//root->addChild(painter->getShapeSep(mesh));
	/*for(int i = 0; i < mesh->verts.size(); i++)
	{
		root->addChild(painter->getSphereSep(mesh, i, 2));
	}*/

	float * color = new float[3];
	color[0] = 0.7f;
	color[1] = 0.7f;
	color[2] = 0.7f;

	float * blue = new float[3];
	blue[0] = 0.0f;
	blue[1] = 0.0f;
	blue[2] = 1.0f;

	if(taskId == 1 || taskId == 2 || taskId == 3)
	{
		/*****************Disk Parameterization******************/
		WeightEnum paramType;
		if(taskId == 1)
		{
			cout << "1: Disk Parameterization | Uniform Weights" << endl;
			paramType = UNIFORM;
		}
		else if(taskId == 2)
		{
			cout << "2: Disk Parameterization | Harmonic Weights" << endl;
			paramType = HARMONIC;
		}
		else if(taskId == 3)
		{
			cout << "3: Disk Parameterization | Mean Value Weights" << endl;
			paramType = MEAN_VALUE;
		}
		
		Parameterization *param = new Parameterization(paramType);
		set<int> boundaryVertices = param->FindBoundaryVertices(mesh);

		vector<FreeVertex *> circleVertices = param->GetCircleVertices(boundaryVertices.size());
		for(int i = 0; i < circleVertices.size() - 1; i++)
		{
			root->addChild(painter->getSphereSepByCoord(mesh, circleVertices[i]->coords, 0.2f));
		}
		root->addChild(painter->getSphereSepByCoord(mesh, circleVertices[circleVertices.size() - 1]->coords, 0.2f));

		vector<FreeVertex *> resultVertices = param->DiskParameterization(mesh);

		// Find the corresponding new vertices for the original vertices wrt their ordering in mesh->edges.
		for(int i = 0; i < mesh->edges.size(); i++)
		{
			Edge *edge = mesh->edges[i];
			float *v1coords = resultVertices[edge->v1i]->coords;
			float *v2coords = resultVertices[edge->v2i]->coords;
			root->addChild(painter->drawEdge(mesh, v1coords, v2coords, color, false));
		}
	}
	else if(taskId == 4)
	{
		/*****************Spherical Parameterization******************/
		cout << "4: Spherical Parameterization of Closed Mesh" << endl;

		bool isConvex = false;
		if(filename[0] == '3' || filename[0] == '7')
		{
			isConvex = true;
		}

		SphereParam *sphereParam = new SphereParam(isConvex);
		FreeVertex *centroid = sphereParam->FindCentroid(mesh);
		/*root->addChild(painter->getSphereSepByCoord(mesh, centroid->coords, 1.0f));
		root->addChild(painter->getShapeSep(mesh));*/
		sphereParam->SphericalParameterization(mesh);
		root->addChild(painter->getShapeSep(mesh));
		for(int i = 0; i < mesh->edges.size(); i++)
		{
			int v1 = mesh->edges[i]->v1i;
			int v2 = mesh->edges[i]->v2i;
			float *v1coords = mesh->verts[v1]->coords;
			float *v2coords = mesh->verts[v2]->coords;
			root->addChild(painter->drawEdge(mesh, v1coords, v2coords, blue, true));
		}
	}
	else if(taskId == 5)
	{
		cout << "5: Disk Parameterization of Closed Mesh" << endl;
		cout << "Not implemented :(" << endl;
	}
	else if(taskId == 6)
	{
		/********************Sphere Generation********************/
		cout << "6: Sphere Generation" << endl;
		
		SphereGen *sphere = new SphereGen();
		Mesh *tetrahedronSphere = sphere->GenerateSphere();

		for(int i = 0; i < tetrahedronSphere->edges.size(); i++)
		{
			Edge *edge = tetrahedronSphere->edges[i];
			float *v1coords = tetrahedronSphere->verts[edge->v1i]->coords;
			float *v2coords = tetrahedronSphere->verts[edge->v2i]->coords;
			root->addChild(painter->drawEdge(tetrahedronSphere, v1coords, v2coords, color, false));
		}
	}
	else if(taskId == 7)
	{
		cout << "7: Disk Parameterization of Generated Sphere" << endl;

		SphereGen *sphere = new SphereGen();
		mesh = sphere->GenerateSphere();
		//mesh->loadOff("sphere.off");

		// Draw the sphere generated from tetrahedron.
		for(int i = 0; i < mesh->edges.size(); i++)
		{
			Edge *edge = mesh->edges[i];
			float *v1coords = mesh->verts[edge->v1i]->coords;
			float *v2coords = mesh->verts[edge->v2i]->coords;
			root->addChild(painter->drawEdge(mesh, v1coords, v2coords, color, false));
		}

		// Draw the pole vertices(cut start and end).
		// vector<int> poleVertices = sphere->FindPoleVertices(mesh);
		// for(int i = 0; i < poleVertices.size(); i++)
		// {
		// 	root->addChild(painter->getSphereSepByCoord(mesh, mesh->verts[poleVertices[i]]->coords, 0.01f));
		// }

		// Draw edges between the cut vertices.
		vector<int> cutVertices;
		sphere->FindCutVertices(mesh, cutVertices);
		for(int i = 0; i < cutVertices.size() - 1; i++)
		{
			float *v1coords = mesh->verts[cutVertices[i]]->coords;
			float *v2coords = mesh->verts[cutVertices[i + 1]]->coords;
			root->addChild(painter->drawEdge(mesh, v1coords, v2coords, blue, true));
		}
		// Draw spheres to cut vertices.
		for(int i = 0; i < cutVertices.size(); i++)
		{
			int v = cutVertices[i];
			root->addChild(painter->getSphereSepByCoord(mesh, mesh->verts[cutVertices[i]]->coords, 0.02f));
			/*vector<int> triList = mesh->verts[v]->triList;
			for(int j  = 0; j < mesh->verts[v]->triList.size(); j++)
			{
				int tri = mesh->verts[v]->triList[j];
				int v1 = mesh->tris[tri]->v1i;
				int v2 = mesh->tris[tri]->v2i;
				int v3 = mesh->tris[tri]->v3i;
				root->addChild(painter->getSphereSep(mesh, v1, 0));
				root->addChild(painter->getSphereSep(mesh, v2, 1));
				root->addChild(painter->getSphereSep(mesh, v3, 2));
			}*/
		}

		//pair<vector<pair<int, int>>, set<int>> cutResult = sphere->CreateCut(mesh, cutVertices);
		//vector<pair<int, int>> triLabels = cutResult.first;
		//set<int> cutTriIds = cutResult.second;
		//cout << cutTriIds.size() << endl;
		//cout << triLabels.size() << endl;
		//cout << cutVertices.size() << endl;

		//for(int i = 0; i < cutTriIds.size(); i++)
		//{
		//	int cutTriId = *next(cutTriIds.begin(), i);
		//	Triangle *tri = mesh->tris[cutTriId];
		//	root->addChild(painter->getSphereSep(mesh, tri->v1i, 1)); //triLabels[i].second));
		//	root->addChild(painter->getSphereSep(mesh, tri->v2i, 1)); //triLabels[i].second));
		//	root->addChild(painter->getSphereSep(mesh, tri->v3i, 1)); //triLabels[i].second));
		//}
		
		// If cut is created successfully
		// set<int> boundaryVertices(cutVertices.begin(), cutVertices.end());
		// Parameterization *p = new Parameterization(UNIFORM);
		// p->ParamClosedMesh(mesh, boundaryVertices);
	}

	// Stuff to be drawn must be added to the root.
	viewer->setSize(SbVec2s(1000, 1000));
	// viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
	return 0;
}
