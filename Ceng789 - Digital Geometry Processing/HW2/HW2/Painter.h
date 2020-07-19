#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoShapeHints.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoIndexedLineSet.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoPointSet.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoTransform.h>
//#include <Inventor/nodes/SoColorMap.h>

#include "Mesh.h"


class Painter
{
public:
	SoSeparator* getShapeSep(Mesh* mesh);
	SoSeparator* drawEdge(Mesh* mesh,  float *v1coords, float *v2coords, float *color, bool isThickEdge);
	SoSeparator* getSphereSep(Mesh* mesh, int v, int i, float radius);
	SoSeparator* getSphereSepByCoord(Mesh* mesh, float *coords, float radius);
	SoSeparator* getSpheresSep(Mesh* mesh, vector<int> samples);
	SoSeparator* getPointSep(Mesh* mesh, float *v1coords);
	SoSeparator* paintVertices(Mesh *mesh, vector<int> samples);
};
