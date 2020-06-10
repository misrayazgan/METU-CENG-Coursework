#include "Sampling.h"

Sampling::Sampling(int n)
{
	numberOfSamples = n;
}

float Sampling::FindLength(float *a)
{
	return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

float* Sampling::CrossProduct(float *a, float *b)
{
	float *result = new float[3];
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];

    return result;
}

float Sampling::CalculateArea(Mesh *mesh, int triId)
{
	// area = 0.5 * |(a-c) x (b-c)|
	Triangle * tri = mesh->tris[triId];
	float *v2tov1 = new float[3];
	float *v2tov3 = new float[3];

	v2tov1[0] = mesh->verts[tri->v1i]->coords[0] - mesh->verts[tri->v2i]->coords[0];
	v2tov1[1] = mesh->verts[tri->v1i]->coords[1] - mesh->verts[tri->v2i]->coords[1];
	v2tov1[2] = mesh->verts[tri->v1i]->coords[2] - mesh->verts[tri->v2i]->coords[2];

	v2tov3[0] = mesh->verts[tri->v3i]->coords[0] - mesh->verts[tri->v2i]->coords[0];
	v2tov3[1] = mesh->verts[tri->v3i]->coords[1] - mesh->verts[tri->v2i]->coords[1];
	v2tov3[2] = mesh->verts[tri->v3i]->coords[2] - mesh->verts[tri->v2i]->coords[2];

	float area = 0.5 * FindLength(CrossProduct(v2tov1, v2tov3));
	return area;
}

vector<float *> Sampling::UniformSampling(Mesh *mesh)
{
	// Pick triangles with probabilities proportional to their areas.
	// Then, put a random sample inside the picked ones.

	// Seed for random number generation
	srand(static_cast<unsigned>(time(0)));

	vector<pair<int, float>> areas;		// (triId, area) pairs
	vector<float> probs(mesh->tris.size(), 0);
	vector<float> cdf(mesh->tris.size(), 0);
	vector<int> selectedTris;
	vector<float *> samples;

	// Find areas of triangles.
	for(int i = 0; i < mesh->tris.size(); i++)
	{
		areas.push_back(make_pair(i, CalculateArea(mesh, i)));
	}

	float sum = accumulate(areas.begin(), areas.end(), 0.0f, [](float &a, pair<int, float> &b){return a + b.second;});

	// Normalize area vector so that it can be used as probabilities.	
	for(int i = 0; i < areas.size(); i++)
	{
		probs[i] = areas[i].second / sum;
	}

	// Calculate CDF (sum the probs one by one (cumulative sum)).
	partial_sum(probs.begin(), probs.end(), cdf.begin());

	for(int n = 0; n < numberOfSamples; n++)
	{
		// Generate a random number between 0 and 1(0.9999.. in our case).
		float random = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX/cdf.back()));

		// Perform binary search to find which interval contains the random number.
		// Use idx as the index to areas vector to obtain picked triangle.
		int idx = std::lower_bound(cdf.begin(), cdf.end(), random) - cdf.begin();
		selectedTris.push_back(areas[idx].first);
	}
	
	vector<int> count(mesh->tris.size(), 0);
	for(int i = 0; i < selectedTris.size(); i++)
	{
		count[selectedTris[i]]++;
	}

	// Now, we have the number of samples per triangle.
	// Sample points inside the triangles by using barycentric coordinates.
	// u + v <= 1, w = 1 - (u+v), S = u * A + v * B + w * C.
	for(int i = 0; i < selectedTris.size(); i++)
	{
		float *a = mesh->verts[mesh->tris[selectedTris[i]]->v1i]->coords;
		float *b = mesh->verts[mesh->tris[selectedTris[i]]->v2i]->coords;
		float *c = mesh->verts[mesh->tris[selectedTris[i]]->v3i]->coords;

		float u, v, w;
		while(true)
		{
			u = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX));
			v = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX));
			if(u + v <= 1)
				break;
		}

		w = 1 - (u + v);

		float *sample = new float[3];
		sample[0] = u * a[0] + v * b[0] + w * c[0];
		sample[1] = u * a[1] + v * b[1] + w * c[1];
		sample[2] = u * a[2] + v * b[2] + w * c[2];

		samples.push_back(sample);
	}
	
	return samples;
}