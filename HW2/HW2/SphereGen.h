#include "Mesh.h"
#include <set>
#include <map>
#include "Dijkstra.h"

class SphereGen
{
public:
	SphereGen() {};
	Mesh * GenerateSphere();
private:
	FreeVertex *centroid;
	FreeVertex* FindCentroid(Mesh *mesh);
	void GetSphereVertices(Mesh *mesh);
	void DivideTriangles(Mesh *mesh, int v1, int v2, int v3, int level);
	float * GetMiddlePoint(float *v1coords, float *v2coords);
	Mesh * CreateTetrahedron();
	void SplitTetrahedron(Mesh *mesh, int level);
};
