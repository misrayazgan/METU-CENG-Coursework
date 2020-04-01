#include <Eigen/Dense>
#include "Mesh.h"
#include <set>
#include <numeric>

#define PI 3.14159265

enum WeightEnum
{
	UNIFORM,
	HARMONIC,
	MEAN_VALUE
};

class Parameterization
{
public:
	Parameterization(WeightEnum w);
	void writeMatrix(Mesh *mesh, Eigen::MatrixXf W);
	vector<FreeVertex *> DiskParameterization(Mesh *mesh);
	set<int> FindBoundaryVertices(Mesh *mesh);
	vector<FreeVertex *> GetCircleVertices(int numberOfVertices);
	vector<FreeVertex *> ClosedMeshParam(Mesh *mesh);
	vector<FreeVertex *> SphereParam(Mesh *mesh, vector<int> boundaryVertices);
private:
	WeightEnum weight;
	vector<bool> isBoundary;
	vector<FreeVertex *> ParamUniform(Mesh *mesh, vector<FreeVertex *> circleVertices);
	vector<FreeVertex *> ParamHarmonic(Mesh *mesh, vector<FreeVertex *> circleVertices);
	vector<FreeVertex *> ParamMeanValue(Mesh *mesh, vector<FreeVertex *> circleVertices);
	int FindClosestCircleVertex(float *srcCoords, vector<FreeVertex *> circleVertices);
	float FindDistance(float* coords1, float* coords2);
	float FindAngle(float *v1coords, float *v2coords);
	float* GetVector(float* p1, float* p2);
	bool isMouthVertex(Mesh *mesh, int v);
	int FindSuccessiveVertex(vector<FreeEdge *> edges, int prevV, int currentV);
	void display(vector<vector<float>> W);
};