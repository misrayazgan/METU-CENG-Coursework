#include <Eigen/Dense>
#include "Mesh.h"
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
	void GetSphereVertices(Mesh *mesh);
	bool isNegated(float *a, float *b);
	float determinant(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2);
	float triangleIntersection(Mesh *mesh, int t, float *vec, float *center);
	Vec3f subtract(const Vec3f &a, const Vec3f &b);
	float * FindIntersectionPoint(float *center, float *vec, float t);
};