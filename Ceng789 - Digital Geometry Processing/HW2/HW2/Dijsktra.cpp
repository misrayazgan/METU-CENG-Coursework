#include "Dijkstra.h"

float Inf = numeric_limits<float>::infinity();

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

	result = DijkstraMinHeap(mesh, sourceIdx);
	vector<int> prev = result.second;

	GetPath(prev, destIdx, path);

	return path;
}
