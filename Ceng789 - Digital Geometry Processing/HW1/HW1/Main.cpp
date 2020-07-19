#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>

#include "Mesh.h"
#include "Painter.h"
#include "Descriptor.h"

int main(int, char ** argv)
{
	HWND window = SoWin::init(argv[0]);

	SoWinExaminerViewer *viewer = new SoWinExaminerViewer(window);

	//make a dead simple scene graph by using the Coin library, only containing a single cone under the scenegraph root
	SoSeparator *root = new SoSeparator;
	root->ref();

	char * filename = "man0.off";
	int taskId;

	cout << "1: Dijkstra | Calculate NxN matrix and write it to file" << endl;
	cout << "2: Dijkstra | Find the shortest path between two query points" << endl;
	cout << "3: Farthest Point Sampling" << endl;
	cout << "4: Gaussian Curvature Descriptor" << endl;
	cout << "5: Average Geodesic Distances Descriptor" << endl;
	cout << "6: Geodesic Isocurve Signature Descriptor" << endl;
	cout << "Enter the id for the desired task:" << endl;
	cin >> taskId;

	DijkstraType type = DIJKSTRA_MINHEAP;

	Mesh *mesh = new Mesh();
	Painter *painter = new Painter();
	Dijkstra *dij = new Dijkstra(type);
	
	mesh->loadOff(filename);
	
	printf("Number of vertices: %d\n", mesh->verts.size());
	
	if(taskId == 1)
	{
		/******** for calculating NxN matrix and writing it to file ********/
		cout << "Dijkstra | Calculate NxN matrix and write it to file" << endl;

		vector<vector<float>> allDist = dij->FindAllDistances(mesh, true);
		root->addChild(painter->getShapeSep(mesh));
	}
	else if(taskId == 2)
	{
		/******** for finding the shortest path between two query points ********/
		cout << "Dijkstra | Find the shortest path between two query points" << endl;

		vector<int> result;
		//auto startTime = chrono::high_resolution_clock::now();
		result = dij->GetShortestPath(mesh, 276, 14938);
		/*auto endTime = chrono::high_resolution_clock::now();
		auto timePassed = endTime - startTime;
		long long microsec = (long long)chrono::duration_cast<chrono::microseconds>(timePassed).count();
		printf("microsec: %ld\n", microsec);*/

		/******** for drawing the shortest path between two query points ********/
		root->addChild(painter->getShapeSep(mesh));
		float *blue = new float[3];
		blue[0] = 0.0f;
		blue[1] = 0.0f;
		blue[2] = 1.0f;
		for(int i = 0; i < result.size() - 1; i++)
		{
			float *v1coords = mesh->verts[result[i]]->coords;
			float *v2coords = mesh->verts[result[i + 1]]->coords;
			SoSeparator * sp = painter->drawThickEdge(mesh, v1coords, v2coords, blue);
			root->addChild(sp);
		}
	}
	else if(taskId == 3)
	{	
		/******** FPS algorithm ********/
		cout << "Farthest Point Sampling" << endl;

		int numberOfSamples = 10;
		Sampling *fps = new Sampling(numberOfSamples);
		vector<int> sampledVertices = fps->FPS(mesh);

		cout << sampledVertices.size() << endl;

		/******** for drawing spheres around sampledVertices ********/
		root->addChild(painter->getShapeSep(mesh));
		root->addChild(painter->getSpheresSep(mesh, sampledVertices));
	}
	else if(taskId == 4)
	{
		/******** Gaussian Curvature Descriptor ********/
		cout << "Gaussian Curvature Descriptor" << endl;

		Descriptor *desc = new Descriptor();
		pair<vector<int>, vector<int>> histogramPair = desc->FindGaussianHistogram(mesh);
		vector<int> histogram = histogramPair.first;
		vector<int> vertsToBinsGaussian = histogramPair.second;
	
		/*for(int i = 0; i < histogram.size(); i++)
		{
			cout << histogram[i] << endl;
		}*/

		root->addChild(painter->paintVertices(mesh, vertsToBinsGaussian));
	}
	else if(taskId == 5)
	{
		/******** Average Geodesic Distances Descriptor ********/
		cout << "Average Geodesic Distances Descriptor" << endl;

		Descriptor *desc = new Descriptor();
		/*vector<vector<float>> allDist = dij->FindAllDistances(mesh, false);
		pair<vector<int>, vector<int>> histogramPair = desc->FindAvgGeoDistHistogram(mesh, allDist);*/
		pair<vector<int>, vector<int>> histogramPair = desc->FindAvgGeoDistHistogramByFPS(mesh);
		vector<int> histogram = histogramPair.first;
		vector<int> vertsToBinsAgd = histogramPair.second;

		/*for(int i = 0; i < histogram.size(); i++)
			cout << histogram[i] << endl;*/

		root->addChild(painter->paintVertices(mesh, vertsToBinsAgd));
	}
	else if (taskId == 6)
	{
		/******** Geodesic Isocurve Signature Descriptor ********/
		cout << "Geodesic Isocurve Signature Descriptor" << endl;

		Descriptor *desc = new Descriptor();
		pair<vector<float>, int> histogramPair = desc->FindIsoCurveHistogram(mesh);
		vector<float> isoCurveHistogram = histogramPair.first;
		int seedVertex = histogramPair.second;

		root->addChild(painter->getShapeSep(mesh));
		float *grey = new float[3];
		grey[0] = 0.2f;
		grey[1] = 0.2f;
		grey[2] = 0.2f;

		vector<vector<pair<Vec3f, Vec3f>>> isoCurvePoints = desc->GetIsoCurvePoints();
		for(int i = 0; i < isoCurvePoints.size(); i++)
		{
			for(int j = 0; j < isoCurvePoints[i].size(); j++)
			{
				Vec3f p1 = isoCurvePoints[i][j].first;
				Vec3f p2 = isoCurvePoints[i][j].second;

				// Convert Vec3f to float*
				float *v1coords = new float[3];
				v1coords[0] = p1.x;
				v1coords[1] = p1.y;
				v1coords[2] = p1.z;

				float *v2coords = new float[3];
				v2coords[0] = p2.x;
				v2coords[1] = p2.y;
				v2coords[2] = p2.z;

				root->addChild(painter->drawThickEdge(mesh, v1coords, v2coords, grey));
			}
		}

		root->addChild(painter->getSphereSep(mesh, seedVertex, 1));
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
