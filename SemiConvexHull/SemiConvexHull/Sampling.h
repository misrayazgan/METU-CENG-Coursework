#pragma once

#include "Mesh.h"
#include <algorithm>
#include <numeric>
#include <functional>
#include <cstdlib>
#include <ctime>

class Sampling
{
private:
	int numberOfSamples;
	float FindLength(float *a);
	float* CrossProduct(float *a, float *b);
	float CalculateArea(Mesh *mesh, int triId);
public:
	Sampling(int n);
	vector<float *> UniformSampling(Mesh *mesh);
};