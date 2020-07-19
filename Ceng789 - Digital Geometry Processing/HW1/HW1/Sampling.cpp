#include "Sampling.h"


Sampling::Sampling(int nSamples)
{
	numberOfSamples = nSamples;
}

vector<int> Sampling::FPS(Mesh *mesh)
{
	vector<int> selectedVertices;
	selectedVertices.push_back(0);

	Dijkstra *d = new Dijkstra(DIJKSTRA_MINHEAP);

	pair<vector<float>, vector<int>> result = d->DijkstraMinHeap(mesh, selectedVertices[0]);
	vector<float> minDistancesFromSelected = result.first;

	for(int i = 0; i < numberOfSamples; i++)
	{
		// Find max distance vertex from first selected vertex(0).
		int maxIdx = max_element(minDistancesFromSelected.begin(), minDistancesFromSelected.end()) - minDistancesFromSelected.begin();
		int maxDistance = *max_element(minDistancesFromSelected.begin(), minDistancesFromSelected.end());

		selectedVertices.push_back(maxIdx);
		pair<vector<float>, vector<int>> res = d->DijkstraMinHeap(mesh, maxIdx);
		vector<float> minDistances = res.first;

		for(int j = 0; j < minDistancesFromSelected.size(); j++)
		{
			if(minDistancesFromSelected[j] > minDistances[j])
			{
				minDistancesFromSelected[j] = minDistances[j];
			}
		}
	}

	selectedVertices.erase(selectedVertices.begin());

	return selectedVertices;
}