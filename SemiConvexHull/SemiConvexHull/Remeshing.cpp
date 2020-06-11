#include "Remeshing.h"

float* Remeshing::FindCrossProduct(float *a, float *b)
{
    float *result = new float[3];
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];

    return result;
}

float Remeshing::FindDotProduct(float *a, float *b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

float Remeshing::FindLength(float *a)
{
	return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

float Remeshing::FindDistance(float *a, float *b)
{
	return sqrt(pow((b[0] - a[0]), 2) + pow((b[1] - a[1]), 2) + pow((b[2] - a[2]), 2));
}

// If v is in list, return its index.
// Else return -1.
int Remeshing::GetIndex(const vector<int> list, int v)
{
	int result = find(list.begin(), list.end(), v) - list.begin();
	if(result < list.size())
	{
		return result; 
	}
	return -1;
}

int Remeshing::FindNN(Mesh *mesh, int v, vector<float *> samples)
{
	int nearestIdx = -1;
	float minDist = numeric_limits<float>::infinity();
	float *src = mesh->verts[v]->coords;

	for(int i = 0; i < samples.size(); i++)
	{
		float dist = FindDistance(src, samples[i]);
		if(dist < minDist)
		{
			minDist = dist;
			nearestIdx = i;
		}
	}

	return nearestIdx;
}

float Remeshing::FindAvgAdjLength(Mesh *mesh, int v)
{
	// Compute average length of the adjacent edges.
	vector<int> edgeList = mesh->verts[v]->edgeList;
	float sum = 0;

	for(int i = 0; i < edgeList.size(); i++)
	{
		sum += mesh->edges[edgeList[i]]->length;
	}
	sum /= edgeList.size();

	return sum;
}

float* Remeshing::FindTriangleCenter(Mesh *mesh, int t)
{
	float *v1 = mesh->verts[mesh->tris[t]->v1i]->coords;
	float *v2 = mesh->verts[mesh->tris[t]->v2i]->coords;
	float *v3 = mesh->verts[mesh->tris[t]->v3i]->coords;
	
	float *center = new float[3];
	center[0] = (v1[0] + v2[0] + v3[0]) / 3;
	center[1] = (v1[1] + v2[1] + v3[1]) / 3;
	center[2] = (v1[2] + v2[2] + v3[2]) / 3;

	return center;
}

float Remeshing::FindMinAngleOfTriangle(float *v1coords, float *v2coords, float *v3coords)
{
	// Find edge lengths.
	float lena = FindDistance(v1coords, v2coords);
	float lenb = FindDistance(v2coords, v3coords);
	float lenc = FindDistance(v1coords, v3coords);

	// Find square of edge lengths.
	float sqlena = pow(lena, 2);
	float sqlenb = pow(lenb, 2);
	float sqlenc = pow(lenc, 2);

	// Cosine law: c^2 = a^2 + b^2 - 2(a)(b)(cos beta)
    float alpha = acos((sqlenb + sqlenc - sqlena) / (2 * lenb * lenc)); 
    float betta = acos((sqlena + sqlenc - sqlenb) / (2 * lena * lenc)); 
    float gamma = acos((sqlena + sqlenb - sqlenc) / (2 * lena * lenb)); 

	// Find min angle
	float minAngle = min(alpha, betta);
	minAngle = min(minAngle, gamma);

	return minAngle;
}

void Remeshing::SubdivideTriangles(Mesh *mesh, float avgLen)
{
	float alpha = sqrt(3);
	bool newTriangleCreated = false;

	for(int k = 0; k < 3; k++)
	{
		vector<Triangle *> tris(mesh->tris.begin(), mesh->tris.end());

		for(int i = 0; i < tris.size(); i++)
		{
			//**********************************************************
			//vector<vector<int>> edgeCount(mesh->verts.size(), vector<int>(mesh->verts.size(), 0));

			//for(int k = 0; k < mesh->tris.size(); k++)
			//{
			//	Triangle *tri = mesh->tris[k];
			//	int i, j;

			//	// 1 and 2
			//	i = min(tri->v1i, tri->v2i);
			//	j = max(tri->v1i, tri->v2i);
			//	edgeCount[i][j]++;

			//	// 2 and 3
			//	i = min(tri->v2i, tri->v3i);
			//	j = max(tri->v2i, tri->v3i);
			//	edgeCount[i][j]++;

			//	// 1 and 3
			//	i = min(tri->v1i, tri->v3i);
			//	j = max(tri->v1i, tri->v3i);
			//	edgeCount[i][j]++;
			//}
			//**********************************************************

			int v1 = tris[i]->v1i;
			int v2 = tris[i]->v2i;
			int v3 = tris[i]->v3i;

			float *center = FindTriangleCenter(mesh, i);
			/*float v1scale = FindAvgAdjLength(mesh, v1);
			float v2scale = FindAvgAdjLength(mesh, v2);
			float v3scale = FindAvgAdjLength(mesh, v3);*/
			float v1scale = avgLen / 3;
			float v2scale = avgLen / 3;
			float v3scale = avgLen / 3;
			float centerScale = (v1scale + v2scale + v3scale) / 3;

			float *v1coords = mesh->verts[v1]->coords;
			float *v2coords = mesh->verts[v2]->coords;
			float *v3coords = mesh->verts[v3]->coords;

			float centerv1Dist = FindDistance(center, v1coords);
			float centerv2Dist = FindDistance(center, v2coords);
			float centerv3Dist = FindDistance(center, v3coords);

			if((alpha * centerv1Dist > centerScale) &&
				(alpha * centerv2Dist > centerScale) &&
				(alpha * centerv3Dist > centerScale) &&
				(alpha * centerv1Dist > v1scale) &&
				(alpha * centerv2Dist > v2scale) &&
				(alpha * centerv3Dist > v3scale))
			{
				/*************Subdivide triangle into 3 little triangles from center**************/
				newTriangleCreated = true;

				// Add center point as a new vertex.
				mesh->addVertex(center[0], center[1], center[2]);
				int centerId = mesh->verts.size() - 1;

				// Add edges for v1-c, v2-c, v3-c.
				mesh->addEdge(v1, centerId);
				mesh->addEdge(v2, centerId);
				mesh->addEdge(v3, centerId);

				// Replace the current big triangle with a little one.
				mesh->tris[i]->v1i = centerId;

				// First triangle is replaced with the big one, it is already adj to v2 and v3.
				// So, just remove big triangle(i) from v1.
				int idx = GetIndex(mesh->verts[v1]->triList, i);
				if(idx > -1)
				{
					mesh->verts[v1]->triList.erase(mesh->verts[v1]->triList.begin() + idx);
				}

				// Add i to triList of center.
				mesh->verts[centerId]->triList.push_back(i);

				// Add other two little triangles.
				int triId1 = mesh->tris.size();
				mesh->addTriangleWithoutEdge(v1, centerId, v3);
				int triId2 = mesh->tris.size();
				mesh->addTriangleWithoutEdge(v1, v2, centerId);

				/*************Relax edges of the big triangle*************/
				float *n = FindPlaneNormal(v1coords, v2coords, v3coords);
				// v1 and v2
				if(/*edgeCount[min(v1,v2)][max(v1, v2)] == 2 &&*/ IsValidEdgeFlip(mesh, v1, v2))
				{
					FlipEdge(mesh, v1, v2, n);
				}

				// v2 and v3
				if(/*edgeCount[min(v2,v3)][max(v2, v3)] == 2 &&*/ IsValidEdgeFlip(mesh, v2, v3))
				{
					FlipEdge(mesh, v2, v3, n);
				}

				// v1 and v3
				if(/*edgeCount[min(v1,v3)][max(v1, v3)] == 2 &&*/ IsValidEdgeFlip(mesh, v1, v3))
				{
					FlipEdge(mesh, v1, v3, n);
				}
			}

		}
	}
}

float* Remeshing::FindPlaneNormal(float *a, float *b, float *c)
{
	// Find vectors ab and ac.
	float *ab = new float[3];
	ab[0] = b[0] - a[0];
	ab[1] = b[1] - a[1];
	ab[2] = b[2] - a[2];

	float *ac = new float[3];
	ac[0] = c[0] - a[0];
	ac[1] = c[1] - a[1];
	ac[2] = c[2] - a[2];

	// Unit normal
	float *planeNormal = FindCrossProduct(ab, ac);
	float len = FindLength(planeNormal);
	planeNormal[0] /= len;
	planeNormal[1] /= len;
	planeNormal[2] /= len;

	return planeNormal;
}

bool Remeshing::IsConvexQuad(Mesh *mesh, int otherVertex, int v1, int v2, int v3)
{
	float *v1coords = mesh->verts[v1]->coords;
	float *v2coords = mesh->verts[v2]->coords;
	float *v3coords = mesh->verts[v3]->coords;
	/*float *v1coords = new float[3];
	v1coords[0] = 3;
	v1coords[1] = 0;
	v1coords[2] = 0;

	float *v2coords = new float[3];
	v2coords[0] = 0;
	v2coords[1] = 0;
	v2coords[2] = 0;

	float *v3coords = new float[3];
	v3coords[0] = 0;
	v3coords[1] = 3;
	v3coords[2] = 0;*/

	float *otherCoords = mesh->verts[otherVertex]->coords;

	// Project otherVertex onto the plane created by v1-v2-v3 triangle.
	float *planeNormal = FindPlaneNormal(v1coords, v2coords, v3coords);

	float *v1toother = new float[3];
	v1toother[0] = otherCoords[0] - v1coords[0];
	v1toother[1] = otherCoords[1] - v1coords[1];
	v1toother[2] = otherCoords[2] - v1coords[2];

	// Dot product gives the distance from point to plane.
	float dist = FindDotProduct(v1toother, planeNormal);

	float *projected = new float[3];
	projected[0] = otherCoords[0] - dist * planeNormal[0];
	projected[1] = otherCoords[1] - dist * planeNormal[1];
	projected[2] = otherCoords[2] - dist * planeNormal[2];

	//***************************************
	/*projected[0] = 2;
	projected[1] = 0;
	projected[2] = 0;*/

	// Find if projected point is inside the triangle.
	float *v1v2 = new float[3];
	v1v2[0] = v2coords[0] - v1coords[0];
	v1v2[1] = v2coords[1] - v1coords[1];
	v1v2[2] = v2coords[2] - v1coords[2];
	
	float *v1v3 = new float[3];
	v1v3[0] = v3coords[0] - v1coords[0];
	v1v3[1] = v3coords[1] - v1coords[1];
	v1v3[2] = v3coords[2] - v1coords[2];

	float *v1projected = new float[3];
	v1projected[0] = projected[0] - v1coords[0];
	v1projected[1] = projected[1] - v1coords[1];
	v1projected[2] = projected[2] - v1coords[2];

	Eigen::MatrixXf A(3, 2);
	Eigen::VectorXf b(3);

	A(0, 0) = v1v2[0];
	A(1, 0) = v1v2[1];
	A(2, 0) = v1v2[2];
	A(0, 1) = v1v3[0];
	A(1, 1) = v1v3[1];
	A(2, 1) = v1v3[2];

	b(0) = v1projected[0];
	b(1) = v1projected[1];
	b(2) = v1projected[2];

	Eigen::ColPivHouseholderQR<Eigen::MatrixXf> dec(A);
	Eigen::VectorXf x = dec.solve(b);

	float mn = x(0);
	float mm = x(1);

	if(x(0) > 0 && x(1) > 0 && (x(0) + x(1) <= 1))
	{
		// Projection of the 4th point is inside the triangle.
		// Not a convex quad.
		return false;
	}

	return true;
}

bool Remeshing::IsValidEdgeFlip(Mesh *mesh, int v1, int v2)
{
	// Find two adjacent triangles to the edge with endpoints v1 and v2.
	vector<int> commonTris = mesh->findCommonTriangles(v1, v2);

	// Find the 4 vertices of the quad formed by two adjacent triangles.
	int vert3 = mesh->getThirdVertexOfTriangle(commonTris[0], v1, v2);
	int vert4 = mesh->getThirdVertexOfTriangle(commonTris[1], v1, v2);

	// Check if quad formed by v1-v2-vert3-vert4 is convex.
	// If not convex, edge flip is illegal.
	bool convex1 = IsConvexQuad(mesh, v1, v2, vert3, vert4);
	bool convex2 = IsConvexQuad(mesh, v2, v1, vert3, vert4);
	bool convex3 = IsConvexQuad(mesh, vert3, v1, v2, vert4);
	bool convex4 = IsConvexQuad(mesh, vert4, v1, v2, vert3);

	if(!(convex1 && convex2 && convex3 && convex4))
	{
		// Quad is not convex.
		return false;
	}

	// Find min angle for triangles v1-v2-vert3, v1-v2-vert4
	float angle1 = FindMinAngleOfTriangle(mesh->verts[v1]->coords, mesh->verts[v2]->coords, mesh->verts[vert3]->coords);
	float angle2 = FindMinAngleOfTriangle(mesh->verts[v1]->coords, mesh->verts[v2]->coords, mesh->verts[vert4]->coords);
	float prevMinAngle = min(angle1, angle2);

	// Find min angle for triangles v1-vert3-vert4, v2-vert3-vert4
	angle1 = FindMinAngleOfTriangle(mesh->verts[v1]->coords, mesh->verts[vert3]->coords, mesh->verts[vert4]->coords);
	angle2 = FindMinAngleOfTriangle(mesh->verts[v2]->coords, mesh->verts[vert3]->coords, mesh->verts[vert4]->coords);
	float currentMinAngle = min(angle1, angle2);

	if(currentMinAngle > prevMinAngle)
	{
		// Flip edge
		return true;				
	}
	return false;
}

void Remeshing::RemoveTriangleInfoFromVertices(Mesh *mesh, int id)
{
	// Remove triangle info from triangle's vertices.
	int v1 = mesh->tris[id]->v1i;
	int v2 = mesh->tris[id]->v2i;
	int v3 = mesh->tris[id]->v3i;

	int idx = GetIndex(mesh->verts[v1]->triList, id);
	if(idx > -1)
	{
		mesh->verts[v1]->triList.erase(mesh->verts[v1]->triList.begin() + idx);
	}

	idx = GetIndex(mesh->verts[v2]->triList, id);
	if(idx > -1)
	{
		mesh->verts[v2]->triList.erase(mesh->verts[v2]->triList.begin() + idx);
	}

	idx = GetIndex(mesh->verts[v3]->triList, id);
	if(idx > -1)
	{
		mesh->verts[v3]->triList.erase(mesh->verts[v3]->triList.begin() + idx);
	}
}

void Remeshing::AddTriangleInfoToVertices(Mesh *mesh, int id)
{
	// Add triangle info to triangle's vertices.
	int v1 = mesh->tris[id]->v1i;
	int v2 = mesh->tris[id]->v2i;
	int v3 = mesh->tris[id]->v3i;

	int idx = GetIndex(mesh->verts[v1]->triList, id);
	if(idx == -1)
	{
		mesh->verts[v1]->triList.push_back(id);
	}

	idx = GetIndex(mesh->verts[v2]->triList, id);
	if(idx == -1)
	{
		mesh->verts[v2]->triList.push_back(id);
	}

	idx = GetIndex(mesh->verts[v3]->triList, id);
	if(idx == -1)
	{
		mesh->verts[v3]->triList.push_back(id);
	}
}

void Remeshing::RemoveEdgeInfoFromVertices(Mesh *mesh, int edgeId)
{
	// Remove edgeId from edgeList of edge vertices.
	int v1 = mesh->edges[edgeId]->v1i;
	int v2 = mesh->edges[edgeId]->v2i;

	int idx = GetIndex(mesh->verts[v1]->edgeList, edgeId);
	if(idx > -1)
	{
		mesh->verts[v1]->edgeList.erase(mesh->verts[v1]->edgeList.begin() + idx);
	}

	idx = GetIndex(mesh->verts[v2]->edgeList, edgeId);
	if(idx > -1)
	{
		mesh->verts[v2]->edgeList.erase(mesh->verts[v2]->edgeList.begin() + idx);
	}
}

void Remeshing::AddEdgeInfoToVertices(Mesh *mesh, int edgeId)
{
	// Add edgeId to edgeList of edge vertices.
	int v1 = mesh->edges[edgeId]->v1i;
	int v2 = mesh->edges[edgeId]->v2i;

	int idx = GetIndex(mesh->verts[v1]->edgeList, edgeId);
	if(idx == -1)
	{
		mesh->verts[v1]->edgeList.push_back(edgeId);
	}

	idx = GetIndex(mesh->verts[v2]->edgeList, edgeId);
	if(idx == -1)
	{
		mesh->verts[v2]->edgeList.push_back(edgeId);
	}
}

void Remeshing::FlipEdge(Mesh *mesh, int v1, int v2, float* normal)
{
	// Find two adjacent triangles to the edge with endpoints v1 and v2.
	vector<int> commonTris = mesh->findCommonTriangles(v1, v2);

	// Find the 4 vertices of the quad formed by two adjacent triangles.
	int vert3 = mesh->getThirdVertexOfTriangle(commonTris[0], v1, v2);
	int vert4 = mesh->getThirdVertexOfTriangle(commonTris[1], v1, v2);

	// Remove triangle info from its vertices' vertLists.
	RemoveTriangleInfoFromVertices(mesh, commonTris[0]);
	RemoveTriangleInfoFromVertices(mesh, commonTris[1]);

	// Replace old triangle vertices with new ones.
	mesh->tris[commonTris[0]]->v1i = v1;
	mesh->tris[commonTris[0]]->v2i = vert3;
	mesh->tris[commonTris[0]]->v3i = vert4;

	float *n = FindPlaneNormal(mesh->verts[v1]->coords, mesh->verts[vert3]->coords, mesh->verts[vert4]->coords);
	if(FindDotProduct(n, normal) <= 0)
	{
		mesh->tris[commonTris[0]]->v1i = vert4;
		mesh->tris[commonTris[0]]->v2i = vert3;
		mesh->tris[commonTris[0]]->v3i = v1;
	}

	mesh->tris[commonTris[1]]->v1i = vert4;
	mesh->tris[commonTris[1]]->v2i = vert3;
	mesh->tris[commonTris[1]]->v3i = v2;

	n = FindPlaneNormal(mesh->verts[vert4]->coords, mesh->verts[vert3]->coords, mesh->verts[v2]->coords);
	if(FindDotProduct(n, normal) <= 0)
	{
		mesh->tris[commonTris[1]]->v1i = v2;
		mesh->tris[commonTris[1]]->v2i = vert3;
		mesh->tris[commonTris[1]]->v3i = vert4;
	}


	// Add triangle info to new triangles' vertices' vertLists.
	AddTriangleInfoToVertices(mesh, commonTris[0]);
	AddTriangleInfoToVertices(mesh, commonTris[1]);

	// Delete v1 and v2 from each other's vertLists.
	int idx = GetIndex(mesh->verts[v1]->vertList, v2);
	if(idx > -1)
	{
		mesh->verts[v1]->vertList.erase(mesh->verts[v1]->vertList.begin() + idx);
	}

	idx = GetIndex(mesh->verts[v2]->vertList, v1);
	if(idx > -1)
	{
		mesh->verts[v2]->vertList.erase(mesh->verts[v2]->vertList.begin() + idx);
	}

	// Remove edge v1-v2 from edgeList of v1 and v2.
	int edgeId = mesh->findEdgeByVertices(v1, v2);
	RemoveEdgeInfoFromVertices(mesh, edgeId);

	// Replace edge v1-v2 with edge vert3-vert4
	mesh->edges[edgeId]->v1i = vert3;
	mesh->edges[edgeId]->v2i = vert4;

	// Add edge vert3-vert4 to edgeList of vert3 and vert4.
	AddEdgeInfoToVertices(mesh, edgeId);
}

// Calculate dE for each vertex.
vector<vector<float>> Remeshing::CalculateEnergyDerivative(Mesh *mesh, vector<float *> samples)
{
	vector<vector<float>> dE(mesh->verts.size(), vector<float>(3, 0));

	// Vertex energy
	for(int i = 0; i < mesh->verts.size(); i++)
	{
		int nnIdx = FindNN(mesh, i, samples);
		
		dE[i][0] = 2 * (mesh->verts[i]->coords[0] - samples[nnIdx][0]);
		dE[i][1] = 2 * (mesh->verts[i]->coords[1] - samples[nnIdx][1]);
		dE[i][2] = 2 * (mesh->verts[i]->coords[2] - samples[nnIdx][2]);
	}

	// Edge energy
	float avgEdgeLen = mesh->computeAvgEdgeLength();
	for(int i = 0; i < mesh->edges.size(); i++)
	{
		int v1 = mesh->edges[i]->v1i;
		int v2 = mesh->edges[i]->v2i;
		float *v1coords = mesh->verts[v1]->coords;
		float *v2coords = mesh->verts[v2]->coords;
		float edgeLen = FindDistance(v1coords, v2coords);
		dE[v1][0] += 4 * (edgeLen * edgeLen - avgEdgeLen) * (v1coords[0] - v2coords[0]);
		dE[v1][1] += 4 * (edgeLen * edgeLen - avgEdgeLen) * (v1coords[1] - v2coords[1]);
		dE[v1][2] += 4 * (edgeLen * edgeLen - avgEdgeLen) * (v1coords[2] - v2coords[2]);

		// dE[v2] lere ne olacak ?????????*
		// dE(v2,:) = dE(v2,:) - 4*(norm2(V(v1,:)-V(v2,:))-l2)*(V(v1,:)-V(v2,:));
	}

	return dE;
}

void Remeshing::ShiftVertices(Mesh *mesh, float gamma, vector<vector<float>> dE)
{
	for(int i = 0; i < mesh->verts.size(); i++)
	{
		mesh->verts[i]->coords[0] -= gamma * dE[i][0];
		mesh->verts[i]->coords[1] -= gamma * dE[i][1];
		mesh->verts[i]->coords[2] -= gamma * dE[i][2];
	}
}

void Remeshing::ScaleVertices(Mesh *mesh, float scale)
{
	for(int i = 0; i < mesh->verts.size(); i++)
	{
		mesh->verts[i]->coords[0] *= scale;
		mesh->verts[i]->coords[1] *= scale;
		mesh->verts[i]->coords[2] *= scale;
	}
}

void Remeshing::OptimizeMesh(Mesh *mesh)
{
	float w = 5;
	float gamma = 0.03;
	float *s = new float[2];
	s[0] = 1.2;
	s[1] = 1.05;
	
	int nSamples = mesh->tris.size() * 2;
	Sampling *sampling = new Sampling(nSamples);
	vector<float* > samples = sampling->UniformSampling(mesh);

	for(int i = 0; i < 15; i++)
	{
		vector<vector<float>> dE = CalculateEnergyDerivative(mesh, samples);
		ShiftVertices(mesh, gamma, dE);
	}

}



