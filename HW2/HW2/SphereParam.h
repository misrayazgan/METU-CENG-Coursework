#include <Eigen/Dense>
#include "Mesh.h"
#include "Sampling.h"
#include <set>
#include <numeric>

class SphereParam
{
public:
	SphereParam(bool convex) { isConvex = convex; };
	void SphericalParameterization(Mesh *mesh);
	FreeVertex * FindCentroid(Mesh *mesh);
private:
	bool isConvex;
	FreeVertex *centroid;
	int FindCollisions(Mesh *mesh);
	//vector<FreeVertex *> RelocateVertices(Mesh * mesh);
	void RelocateVertices(Mesh * mesh);
	float FindDistance(float* coords1, float* coords2);
	float FindVectorLen(float *vec);
	void GetSphereVertices(Mesh *mesh);
	bool isNegated(float *a, float *b);
	float Determinant(float *v0, float *v1, float *v2);
	float TriangleIntersection(Mesh *mesh, int t, float *vec, float *center);
	float * Subtract(float *a, float *b);
	float * FindIntersectionPoint(float *center, float *vec, float t);
};
