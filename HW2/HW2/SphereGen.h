#include "Mesh.h"
#include <set>
#include "Dijkstra.h"

class SphereGen
{
public:
	SphereGen() {};
	Mesh * GenerateSphere();
	vector<int> FindPoleVertices(Mesh *mesh);
	vector<int> FindCutVertices(Mesh *mesh);
	pair<vector<int>, set<int>> CreateCut(Mesh *mesh);
private:
	FreeVertex *centroid;
	FreeVertex* FindCentroid(Mesh *mesh);
	void GetSphereVertices(Mesh *mesh);
	void DivideTriangles(Mesh *mesh, int v1, int v2, int v3, int level);
	float * GetMiddlePoint(float *v1coords, float *v2coords);
	Mesh * CreateTetrahedron();
	void SplitTetrahedron(Mesh *mesh, int level);
	
	vector<int> FindCommonVertices(Mesh *mesh, int t1, int t2);
	bool IsCutEdge(vector<int> cutVertices, vector<int> verts);
};