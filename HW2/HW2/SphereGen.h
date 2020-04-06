#include "Mesh.h"
#include <set>
#include <map>
#include "Dijkstra.h"

class SphereGen
{
public:
	SphereGen() {};
	Mesh * GenerateSphere();
	vector<int> FindPoleVertices(Mesh *mesh);
	void FindCutVertices(Mesh *mesh, vector<int> &cutVertices);
	pair<map<int, int>, set<int>> CreateCut(Mesh *mesh, vector<int> &cutVertices, vector<int> &duplicateVertices);
	int GetIndex(const vector<int> cutVertices, int v);
private:
	FreeVertex *centroid;
	FreeVertex* FindCentroid(Mesh *mesh);
	void GetSphereVertices(Mesh *mesh);
	void DivideTriangles(Mesh *mesh, int v1, int v2, int v3, int level);
	float * GetMiddlePoint(float *v1coords, float *v2coords);
	Mesh * CreateTetrahedron();
	void SplitTetrahedron(Mesh *mesh, int level);
	void FillLabels(map<int, int> triLabels, vector<int> &labeled0, vector<int> &labeled1);
	void FillCutsNonCuts(Mesh *mesh, vector<int> &cuts, vector<int> &noncuts, vector<int> cutVertices, int triId);
	
	vector<int> FindCommonVertices(Mesh *mesh, int t1, int t2);
	bool IsCutEdge(vector<int> cutVertices, vector<int> verts);
};
