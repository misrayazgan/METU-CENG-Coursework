#include "Descriptor.h"

vector<float> Descriptor::FindGaussianCurvature(Mesh *mesh)
{
	vector<float> allK;

	// Iterate over all vertices of the mesh.
	vector<Vertex *>::iterator it;
	for(it = mesh->verts.begin(); it != mesh->verts.end(); ++it)
	{
		Vertex *v = *it;
		float k = 360.0;

		// Consider all adjacent triangles to the vertex.
		for(int i = 0; i < v->triList.size(); i++)
		{
			float theta = FindAngle(mesh, v, v->triList[i]);
			k -= theta;
		}

		allK.push_back(k);
	}

	return allK;
}

pair<vector<int>, vector<int>> Descriptor::FindGaussianHistogram(Mesh *mesh)
{
	vector<float> allK = FindGaussianCurvature(mesh);

	float min = *min_element(allK.begin(), allK.end());
	float max = *max_element(allK.begin(), allK.end());

	// Add 1 to include max
	//float range = max - min + 1;
	//int numberOfBins = 5;
	//float stepSize = range / numberOfBins;

	//// Store number of vertices in histogram and the corresponding vertices in another vector.
	//vector<int> hist(numberOfBins, 0);
	//vector<int> vertsToBins(mesh->verts.size(), 0);

	//for(int v = 0; v < allK.size(); v++)
	//{
	//	float k = allK[v];
	//	int bin = (int)((k - min)/stepSize);
	//	hist[bin]++;
	//	vertsToBins[v] = bin;
	//}

	int numberOfBins = 3;
	vector<int> hist(numberOfBins, 0);
	vector<int> vertsToBins(mesh->verts.size(), 0);

	for(int v = 0; v < allK.size(); v++)
	{
		float k = allK[v];
		int bin;
		if(k < 0)
			bin = 0;
		else if(k == 0)
			bin = 1;
		else
			bin = 2;

		hist[bin]++;
		vertsToBins[v] = bin;
	}

	return make_pair(hist, vertsToBins);
}

float Descriptor::FindAngle(Mesh *mesh, Vertex *v, int triId)
{
	// x.y = |x|.|y|.cosa
	float theta;
	Triangle *tri = mesh->tris[triId];
	int v1i = tri->v1i;
	int v2i = tri->v2i;
	int v3i = tri->v3i;

	Vertex *v1 = mesh->verts[v1i];
	Vertex *v2 = mesh->verts[v2i];
	Vertex *v3 = mesh->verts[v3i];

	Vec3f xVec;
	Vec3f yVec;

	if(v->idx == v1i)
	{
		float xLen = mesh->computeEdgeLength(v1i, v2i);
		float yLen = mesh->computeEdgeLength(v1i, v3i);	
		
		xVec.x = v2->coords[0] - v1->coords[0];
		xVec.y = v2->coords[1] - v1->coords[1];
		xVec.z = v2->coords[2] - v1->coords[2];

		yVec.x = v3->coords[0] - v1->coords[0];
		yVec.y = v3->coords[1] - v1->coords[1];
		yVec.z = v3->coords[2] - v1->coords[2];
		
		float dotPro = ComputeDotProduct(xVec, yVec);
		theta = acos(dotPro / (xLen * yLen)) * 180.0 / PI;
	}
	else if(v->idx == v2i)
	{
		float xLen = mesh->computeEdgeLength(v2i, v1i);
		float yLen = mesh->computeEdgeLength(v2i, v3i);	
		
		xVec.x = v1->coords[0] - v2->coords[0];
		xVec.y = v1->coords[1] - v2->coords[1];
		xVec.z = v1->coords[2] - v2->coords[2];

		yVec.x = v3->coords[0] - v2->coords[0];
		yVec.y = v3->coords[1] - v2->coords[1];
		yVec.z = v3->coords[2] - v2->coords[2];
		
		float dotPro = ComputeDotProduct(xVec, yVec);
		theta = acos(dotPro / (xLen * yLen)) * 180.0 / PI;
	}
	else if(v->idx == v3i)
	{
		float xLen = mesh->computeEdgeLength(v3i, v2i);
		float yLen = mesh->computeEdgeLength(v3i, v1i);	
		
		xVec.x = v2->coords[0] - v3->coords[0];
		xVec.y = v2->coords[1] - v3->coords[1];
		xVec.z = v2->coords[2] - v3->coords[2];

		yVec.x = v1->coords[0] - v3->coords[0];
		yVec.y = v1->coords[1] - v3->coords[1];
		yVec.z = v1->coords[2] - v3->coords[2];
		
		float dotPro = ComputeDotProduct(xVec, yVec);
		theta = acos(dotPro / (xLen * yLen)) * 180.0 / PI;
	}

	return theta;
}

float Descriptor::ComputeDotProduct(Vec3f xVec, Vec3f yVec)
{
	float result = xVec.x * yVec.x + xVec.y * yVec.y + xVec.z * yVec.z;
	return result;
}

vector<float> Descriptor::FindAverageGeodesicDistances(Mesh *mesh, vector<vector<float>> allDist)
{
	vector<float> result;

	for(int v = 0; v < mesh->verts.size(); v++)
	{
		// Find sum of geodesics from v to all other vertices.
		vector<float> fromV = allDist[v];
		float sum = accumulate(fromV.begin(), fromV.end(), 0);
		result.push_back(sum);
	}
    
    return result;
}

vector<float> Descriptor::FindAverageGeodesicDistancesByFPS(Mesh *mesh)
{
	vector<float> result;

	// Sample 100 vertices
	int numberOfSamples = 100;
	Sampling *fps = new Sampling(numberOfSamples);
	vector<int> sampledVertices = fps->FPS(mesh);

	Dijkstra *dij = new Dijkstra(DIJKSTRA_MINHEAP);
	vector<vector<float>> fromSamplesToOthers;

	for(int i = 0; i < numberOfSamples; i++)
	{
		pair<vector<float>, vector<int>> dijResult = dij->DijkstraMinHeap(mesh, sampledVertices[i]);
		fromSamplesToOthers.push_back(dijResult.first);
	}

	// AGD of a non-sample vertex v is found using the geodesic distances from the samples to v.
	for(int v = 0; v < mesh->verts.size(); v++)
	{
		// Find sum of geodesics from samples to v.
		float sumDist = 0;
		for(int j = 0; j < numberOfSamples; j++)
		{
			sumDist += fromSamplesToOthers[j][v];
		}
		result.push_back(sumDist);
	}

    return result;
}

pair<vector<int>, vector<int>> Descriptor::FindAvgGeoDistHistogramByFPS(Mesh *mesh)
{
	vector<float> sumDist = FindAverageGeodesicDistancesByFPS(mesh);

	float min = *min_element(sumDist.begin(), sumDist.end());
	float max = *max_element(sumDist.begin(), sumDist.end());

	// Add 1 to include max
	float range = max - min + 1;
	int numberOfBins = 5;
	float stepSize = range / numberOfBins;

	// Store number of vertices in histogram and the corresponding vertices in another vector.
	vector<int> hist(numberOfBins, 0);
	vector<int> vertsToBins(mesh->verts.size(), 0);

	for(int v = 0; v < mesh->verts.size(); v++)
	{
		int bin = (int)((sumDist[v] - min) / stepSize);
		hist[bin]++;
		vertsToBins[v] = bin;
	}

	return make_pair(hist, vertsToBins);
}

pair<vector<int>, vector<int>> Descriptor::FindAvgGeoDistHistogram(Mesh *mesh, vector<vector<float>> allDist)
{
	vector<float> sumDist = FindAverageGeodesicDistances(mesh, allDist);

	float min = *min_element(sumDist.begin(), sumDist.end());
	float max = *max_element(sumDist.begin(), sumDist.end());

	// Add 1 to include max
	float range = max - min + 1;
	int numberOfBins = 5;
	float stepSize = range / numberOfBins;

	// Store number of vertices in histogram and the corresponding vertices in another vector.
	vector<int> hist(numberOfBins, 0);
	vector<int> vertsToBins(mesh->verts.size(), 0);

	for(int v = 0; v < mesh->verts.size(); v++)
	{
		int bin = (int)((sumDist[v] - min) / stepSize);
		hist[bin]++;
		vertsToBins[v] = bin;
	}

	return make_pair(hist, vertsToBins);
}

pair<vector<float>, int> Descriptor::FindIsoCurveHistogram(Mesh *mesh)
{
	// Store lengths of isocurves for k different radii (k = numberOfIsocurves)
	vector<float> isoCurveHistogram;

	// Select a seed vertex
	Sampling *fps = new Sampling(10);
	vector<int> sample = fps->FPS(mesh);
	int seed = sample[0];
	//int seed = 0;

	// Find the max radius for isocurves
	Dijkstra *dij = new Dijkstra(DIJKSTRA_MINHEAP);
	pair<vector<float>, vector<int>> dijResult = dij->DijkstraMinHeap(mesh, seed);
	vector<float> distFromSeed = dijResult.first;
	float maxRadius = *max_element(distFromSeed.begin(), distFromSeed.end());
	int numberOfIsocurves = 20;

	// Find radii from r1 ro rn
	float minRadius = maxRadius / numberOfIsocurves;
	vector<float> radii;
	for(int i = 1; i <= numberOfIsocurves; i++)
	{
		radii.push_back(minRadius * i);
	}

	for(int c = 0; c < numberOfIsocurves; c++)
	{
		float radius = radii[c];
		float isoCurveLen = 0;
		vector<pair<Vec3f, Vec3f>> currentCurvePoints;

		// For each triangle
		for(int i = 0; i < mesh->tris.size(); i++)
		{
			// Consider all 3 edges for the triangle. If isocurve intersects two edges, find the common vertex and go on wrt that.
			Triangle *tri = mesh->tris[i];
			int v1b = 0, v2b = 0, v3b = 0;

			// Edge between v1 and v2
			if((distFromSeed[tri->v1i] > radius && distFromSeed[tri->v2i] < radius) || 
				(distFromSeed[tri->v1i] < radius && distFromSeed[tri->v2i] > radius))
			{
				v1b++;
				v2b++;
			}

			// Edge between v2 and v3
			if((distFromSeed[tri->v2i] > radius && distFromSeed[tri->v3i] < radius) || 
				(distFromSeed[tri->v2i] < radius && distFromSeed[tri->v3i] > radius))
			{
				v2b++;
				v3b++;
			}

			// Edge between v1 and v3
			if((distFromSeed[tri->v1i] > radius && distFromSeed[tri->v3i] < radius) || 
				(distFromSeed[tri->v1i] < radius && distFromSeed[tri->v3i] > radius))
			{
				v1b++;
				v3b++;
			}

			if(v1b == 2)
			{
				pair<Vec3f, Vec3f> points = FindPointsOnEdges(mesh, tri->v1i, tri->v2i, tri->v3i, radius, distFromSeed);
				isoCurveLen += ComputeDistance(points.first, points.second);
				currentCurvePoints.push_back(points);
			}
			else if(v2b == 2)
			{
				pair<Vec3f, Vec3f> points = FindPointsOnEdges(mesh, tri->v2i, tri->v1i, tri->v3i, radius, distFromSeed);
				isoCurveLen += ComputeDistance(points.first, points.second);
				currentCurvePoints.push_back(points);
			}
			else if(v3b == 2)
			{
				pair<Vec3f, Vec3f> points = FindPointsOnEdges(mesh, tri->v3i, tri->v1i, tri->v2i, radius, distFromSeed);
				isoCurveLen += ComputeDistance(points.first, points.second);
				currentCurvePoints.push_back(points);
			}
		}
		
		isoCurvePoints.push_back(currentCurvePoints);
		isoCurveHistogram.push_back(isoCurveLen);
		cout << c << " iso curve length is: " << isoCurveLen << endl;
	}

	return make_pair(isoCurveHistogram, seed);
}

vector<vector<pair<Vec3f, Vec3f>>> Descriptor::GetIsoCurvePoints()
{
	return isoCurvePoints;
}

// v0 is  the common vertex (for edges v0-v1 and v0-v2)
pair<Vec3f, Vec3f> Descriptor::FindPointsOnEdges(Mesh *mesh, int v0, int v1, int v2, float radius, vector<float> distFromSeed)
{
	float alpha1 = abs(radius - distFromSeed[v0]) / abs(distFromSeed[v1] - distFromSeed[v0]);
	float alpha2 = abs(radius - distFromSeed[v0]) / abs(distFromSeed[v2] - distFromSeed[v0]);

	Vec3f p1;
	p1.x = (1 - alpha1) * mesh->verts[v0]->coords[0] + alpha1 * mesh->verts[v1]->coords[0];
	p1.y = (1 - alpha1) * mesh->verts[v0]->coords[1] + alpha1 * mesh->verts[v1]->coords[1];
	p1.z = (1 - alpha1) * mesh->verts[v0]->coords[2] + alpha1 * mesh->verts[v1]->coords[2];

	Vec3f p2;
	p2.x = (1 - alpha2) * mesh->verts[v0]->coords[0] + alpha2 * mesh->verts[v2]->coords[0];
	p2.y = (1 - alpha2) * mesh->verts[v0]->coords[1] + alpha2 * mesh->verts[v2]->coords[1];
	p2.z = (1 - alpha2) * mesh->verts[v0]->coords[2] + alpha2 * mesh->verts[v2]->coords[2];

	return make_pair(p1, p2);
}

float Descriptor::ComputeDistance(Vec3f a, Vec3f b)
{
	float distance = sqrt(pow((b.x - a.x), 2) + pow((b.y - a.y), 2) + pow((b.z - a.z), 2));
	return distance;
}

