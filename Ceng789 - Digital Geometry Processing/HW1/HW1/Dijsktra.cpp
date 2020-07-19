#include "Dijkstra.h"

float Inf = numeric_limits<float>::infinity();
char *OutputFilename = "man0_matrix.txt";

Dijkstra::Dijkstra(DijkstraType t)
{
	type = t;
	/*if(t == DIJKSTRA_FIBOHEAP)
	{
		numberOfFiboNodes = 0;
		minPtr = NULL;
	}*/
}

// Dijkstra with array
pair<vector<float>, vector<int>> Dijkstra::DijkstraArray(const Mesh *mesh, int sourceIdx)
{
	Vertex *sourceVert = mesh->verts[sourceIdx];
	int numberOfVertices = mesh->verts.size();

	// Store (distance, vertId) pairs that are being processed in a vector.
	vector<pair<float, int>> pqVec;
	// Init all distances to infinity.
	vector<float> minDistances(numberOfVertices, Inf);
	// Store the previous vertexId for each vertex.
	vector<int> prev(numberOfVertices, 0);

	// Insert source to pq.
	pqVec.push_back(make_pair(0, sourceIdx));
	minDistances[sourceIdx] = 0;
	prev[sourceIdx] = -1;

	// Process vertices until priority queue is empty.
	while(!pqVec.empty())
	{
		// Take the vertex with min distance.
		int minDistId = pqVec[0].second;	// index for verts in mesh
		Vertex * minDistVertex = mesh->verts[minDistId];
		
		// Erase the first element (with min distance)
		pqVec.erase(pqVec.begin());

		// Consider all adjacent vertices for minDistVertex.
		vector<int>::iterator it = (minDistVertex->vertList).begin();
		int i = 0;
		for(; it != (minDistVertex->vertList).end(); ++it)
		{
			int vertIdx = *it;	// index for verts in mesh
			int eIdx = i;	// index for edgeList in vertex
			int edgeIdx = minDistVertex->edgeList[i];	// index for edges in mesh

			float edgeLen = mesh->edges[edgeIdx]->length;

			if (minDistances[minDistId] + edgeLen < minDistances[vertIdx])
			{
				// Shorter path found, update the min distances.
				minDistances[vertIdx] = minDistances[minDistId] + edgeLen;
				pqVec.push_back(make_pair(minDistances[vertIdx], vertIdx));
				SortVector(pqVec);
				// Update prev
				prev[vertIdx] = minDistId;
			}

			i++;
		}
	}

	return make_pair(minDistances, prev);
}

void Dijkstra::SortVector(vector<pair<float,int>> &v)
{
	Functor f;

	// Sort entries in increasing order of distances.
	sort(v.begin(), v.end(), f);
}

// Dijkstra with min heap
pair<vector<float>, vector<int>> Dijkstra::DijkstraMinHeap(const Mesh *mesh, int sourceIdx)
{
	Vertex *sourceVert = mesh->verts[sourceIdx];
	int numberOfVertices = mesh->verts.size();

	// Store (distance, vertId) pairs that are being processed in a priority queue.
    priority_queue<pair<float,int>, vector<pair<float,int>>, greater<pair<float,int>>> pq;
	// Init all distances to infinity.
	vector<float> minDistances(numberOfVertices, Inf);
	// Store the previous vertexId for each vertex.
	vector<int> prev(numberOfVertices, 0);

	// Insert source to pq.
	pq.push(make_pair(0, sourceIdx));
	minDistances[sourceIdx] = 0;
	prev[sourceIdx] = -1;

	// Process vertices until priority queue is empty.
	while(!pq.empty())
	{
		// Take the vertex with min distance.
		int minDistId = pq.top().second;	// index for verts in mesh
		Vertex * minDistVertex = mesh->verts[minDistId];
		pq.pop();

		// Consider all adjacent vertices for minDistVertex.
		vector<int>::iterator it = (minDistVertex->vertList).begin();
		int i = 0;
		for(; it != (minDistVertex->vertList).end(); ++it)
		{
			int vertIdx = *it;	// index for verts in mesh
			int eIdx = i;	// index for edgeList in vertex
			int edgeIdx = minDistVertex->edgeList[i];	// index for edges in mesh

			float edgeLen = mesh->edges[edgeIdx]->length;
			//printf("%d %f\n", *it, edgeLen);

			if (minDistances[minDistId] + edgeLen < minDistances[vertIdx])
			{
				// Shorter path found, update the min distances.
				minDistances[vertIdx] = minDistances[minDistId] + edgeLen;
				pq.push(make_pair(minDistances[vertIdx], vertIdx));
				// Update prev
				prev[vertIdx] = minDistId;
			}

			i++;
		}
	}

	return make_pair(minDistances, prev);
}

// Dijkstra with Fibonacci heap
pair<vector<float>, vector<int>> Dijkstra::DijkstraFiboHeap(const Mesh *mesh, int sourceIdx)
{
	Vertex *sourceVert = mesh->verts[sourceIdx];
	int numberOfVertices = mesh->verts.size();

	// Store (distance, vertId) pairs that are being processed in a fibonacci heap.


	// Init all distances to infinity.
	vector<float> minDistances(numberOfVertices, Inf);
	// Store the previous vertexId for each vertex.
	vector<int> prev(numberOfVertices, 0);

 //   // Insert source to fibo heap.
	//FiboHeapInsertNode(sourceIdx, 0);
	//minDistances[sourceIdx] = 0;
	//prev[sourceIdx] = -1;

	//// Process vertices until priority queue is empty.
	//while(minPtr != NULL)
	//{
	//	// Take the vertex with min distance.
	//	pair<int, float> result = FiboHeapDeleteMin();
	//	int minDistId = result.first;	// index for verts in mesh
	//	float minDist = result.second;

	//	Vertex *minDistVertex = mesh->verts[minDistId];
	//	
	//	// Consider all adjacent vertices for minDistVertex.
	//	vector<int>::iterator it = (minDistVertex->vertList).begin();
	//	int i = 0;
	//	for(; it != (minDistVertex->vertList).end(); ++it)
	//	{
	//		int vertIdx = *it;	// index for verts in mesh
	//		int eIdx = i;	// index for edgeList in vertex
	//		int edgeIdx = minDistVertex->edgeList[i];	// index for edges in mesh

	//		float edgeLen = mesh->edges[edgeIdx]->length;
	//		if (minDistances[minDistId] + edgeLen < minDistances[vertIdx])
	//		{
	//			// Shorter path found, update the min distances.
	//			minDistances[vertIdx] = minDistances[minDistId] + edgeLen;

	//			//Node *
	//			//FiboHeapDecreaseKey();
	//			//pq.push(make_pair(minDistances[vertIdx], vertIdx));
	//			// Update prev
	//			prev[vertIdx] = minDistId;
	//		}

	//		i++;
	//	}
	//}

	return make_pair(minDistances, prev);
}

vector<vector<float>> Dijkstra::FindAllDistances(Mesh *mesh, bool write)
{
	int numberOfVertices = mesh->verts.size();
	vector<vector<float>> allDist;

	if (type == DIJKSTRA_ARRAY)
	{
		vector<Vertex*>::iterator it = (mesh->verts).begin();
		for(; it != (mesh->verts).end(); ++it)
		{
			Vertex *source = *it;
			pair<vector<float>, vector<int>> result = DijkstraArray(mesh, source->idx);
			vector<float> minDistances = result.first;
			vector<int> prev = result.second;
			allDist.push_back(minDistances);
		}
	}
	else if(type == DIJKSTRA_MINHEAP)
	{
		vector<Vertex*>::iterator it = (mesh->verts).begin();
		for(; it != (mesh->verts).end(); ++it)
		{
			Vertex *source = *it;
			pair<vector<float>, vector<int>> result = DijkstraMinHeap(mesh, source->idx);
			vector<float> minDistances = result.first;
			vector<int> prev = result.second;
			allDist.push_back(minDistances);
		}
	}
	else if(type == DIJKSTRA_FIBOHEAP)
	{
		/*vector<Vertex*>::iterator it = (mesh->verts).begin();
		for(; it != (mesh->verts).end(); ++it)
		{
			Vertex *source = *it;
			pair<vector<float>, vector<int>> result = DijkstraFiboHeap(mesh, source->idx);
			vector<float> minDistances = result.first;
			vector<int> prev = result.second;
			allDist.push_back(minDistances);
		}*/
	}	

	if(write)
	{
		WriteMatrix(OutputFilename, allDist);
	}

	return allDist;
}

void Dijkstra::WriteMatrix(char* filename, vector<vector<float>> matrix)
{
	FILE* file = fopen(filename, "w");

	if(file != NULL)
	{
		vector<vector<float>>::iterator i;
		for(i = matrix.begin(); i < matrix.end(); ++i)
		{
			vector<float> row = *i;
			vector<float>::iterator j;
			for(j = row.begin(); j < row.end(); ++j)
			{
				fprintf(file, "%.6g ", *j);
			}
			fprintf(file, "\n");
		}

		fclose(file);
	}
}


void Dijkstra::GetPath(vector<int> prev, int v, vector<int> &result)
{
	// If v is the source
	if(prev[v] == -1)
	{
		result.push_back(v);
		return;
	}
	GetPath(prev, prev[v], result);
	result.push_back(v);
}

vector<int> Dijkstra::GetShortestPath(Mesh *mesh, int sourceIdx, int destIdx)
{
	int numberOfVertices = mesh->verts.size();
	vector<int> path;
	pair<vector<float>, vector<int>> result;

	if(type == DIJKSTRA_ARRAY)
	{
		result = DijkstraArray(mesh, sourceIdx);
	}
	else if(type == DIJKSTRA_MINHEAP)
	{
		result = DijkstraMinHeap(mesh, sourceIdx);
	}
	else if(type == DIJKSTRA_FIBOHEAP)
	{
		result = DijkstraFiboHeap(mesh, sourceIdx);
	}

	vector<int> prev = result.second;

	GetPath(prev, destIdx, path);

	return path;
}
