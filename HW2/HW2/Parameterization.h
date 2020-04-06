#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Mesh.h"
#include "Dijkstra.h"
#include <set>
#include <map>
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
	Parameterization(WeightEnum w, bool closed);
	vector<FreeVertex *> DiskParameterization(Mesh *mesh);
	set<int> FindBoundaryVertices(Mesh *mesh);
	vector<FreeVertex *> GetCircleVertices(int numberOfVertices);

	vector<FreeVertex *> ParamClosedMesh(Mesh *mesh, vector<int> boundaryVertices);
	vector<int> FindPoleVertices(Mesh *mesh);
	void FindCutVertices(Mesh *mesh, vector<int> &cutVertices);
	int GetIndex(const vector<int> cutVertices, int v);
	pair<map<int, int>, set<int>> CreateCut(Mesh *mesh, vector<int> &cutVertices, vector<int> &duplicateVertices);
private:
	WeightEnum weight;
	bool isClosed;
	vector<bool> isBoundary;
	vector<FreeVertex *> ParamUniform(Mesh *mesh, vector<FreeVertex *> circleVertices, vector<int> boundaryVerts);
	vector<FreeVertex *> ParamHarmonic(Mesh *mesh, vector<FreeVertex *> circleVertices);
	vector<FreeVertex *> ParamMeanValue(Mesh *mesh, vector<FreeVertex *> circleVertices);
	int FindClosestCircleVertex(float *srcCoords, vector<FreeVertex *> circleVertices);
	float FindDistance(float* coords1, float* coords2);
	float FindAngle(float *v1coords, float *v2coords);
	float* GetVector(float* p1, float* p2);
	int FindSuccessiveVertex(vector<FreeEdge *> edges, int prevV, int currentV);
	void display(vector<vector<float>> W);
	void SetBoundaryVertices(Mesh *mesh, set<int> boundaryVertices);
	bool IsCutEdge(vector<int> cutVertices, vector<int> verts);
	vector<int> FindCommonVertices(Mesh *mesh, int t1, int t2);
	void FillCutsNonCuts(Mesh *mesh, vector<int> &cuts, vector<int> &noncuts, vector<int> cutVertices, int triId);
	void FillLabels(map<int, int> triLabels, vector<int> &labeled0, vector<int> &labeled1);
};
