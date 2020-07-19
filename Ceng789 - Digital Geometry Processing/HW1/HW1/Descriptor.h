#include "Mesh.h"
#include "Sampling.h"
#include <numeric>
#include <algorithm>
#include <math.h>

#define PI 3.14159265

struct Vec3f
{
	float x;
	float y;
	float z;
};

class Descriptor
{
public:
	vector<float> FindGaussianCurvature(Mesh *mesh);
	pair<vector<int>, vector<int>> FindGaussianHistogram(Mesh *mesh);
	vector<float> FindAverageGeodesicDistances(Mesh *mesh, vector<vector<float>> allDist);
	pair<vector<int>, vector<int>> FindAvgGeoDistHistogram(Mesh *mesh, vector<vector<float>> allDist);
	vector<float> FindAverageGeodesicDistancesByFPS(Mesh *mesh);	
	pair<vector<int>, vector<int>> FindAvgGeoDistHistogramByFPS(Mesh *mesh);	
	pair<vector<float>, int> FindIsoCurveHistogram(Mesh *mesh);
	vector<vector<pair<Vec3f, Vec3f>>> GetIsoCurvePoints();
	
private:
	float FindAngle(Mesh *mesh, Vertex *v, int triId);
	float ComputeDotProduct(Vec3f xVec, Vec3f yVec);
	pair<Vec3f, Vec3f> FindPointsOnEdges(Mesh *mesh, int v0, int v1, int v2, float radius, vector<float> distFromSeed);
	float ComputeDistance(Vec3f a, Vec3f b);
	vector<vector<pair<Vec3f, Vec3f>>> isoCurvePoints;
};