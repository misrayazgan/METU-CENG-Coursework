#include "SphereParam.h"

FreeVertex* SphereParam::FindCentroid(Mesh *mesh)
{
	float *coords = new float[3];
	coords[0] = 0;
	coords[1] = 0;
	coords[2] = 0;

	for(int i = 0; i < mesh->verts.size(); i++)
	{
		coords[0] += mesh->verts[i]->coords[0];
		coords[1] += mesh->verts[i]->coords[1];
		coords[2] += mesh->verts[i]->coords[2];
	}

	coords[0] /= mesh->verts.size();
	coords[1] /= mesh->verts.size();
	coords[2] /= mesh->verts.size();

	if(isConvex == true)
		return new FreeVertex(coords);
	else
	{
		// Find the closest vertex to found center.
		float minDist = numeric_limits<float>::infinity();;
		int minV;
		for(int i = 0; i < mesh->verts.size(); i++)
		{
			float dist = FindDistance(coords, mesh->verts[i]->coords);
			if(dist < minDist)
			{
				minDist = dist;
				minV = i;
			}
		}

		float *closestCoords = mesh->verts[minV]->coords;
		float *ray = new float[3];
		ray[0] = closestCoords[0] - coords[0];
		ray[1] = closestCoords[1] - coords[1];
		ray[2] = closestCoords[2] - coords[2];

		float *concaveCenter = new float[3];

		// Triangle intersection
		for(int i = 0; i < mesh->tris.size(); i++)
		{
			if(mesh->tris[i]->v1i != minV && mesh->tris[i]->v2i != minV && mesh->tris[i]->v3i != minV)
			{
				float hitResult = triangleIntersection(mesh, i, ray, coords);
				if(hitResult > 0)
				{
					// Find intersection point on the triangle.
					float *hitPoint = FindIntersectionPoint(coords, ray, hitResult);
					concaveCenter[0] = (closestCoords[0] + hitPoint[0]) / 2;
					concaveCenter[1] = (closestCoords[1] + hitPoint[1]) / 2;
					concaveCenter[2] = (closestCoords[2] + hitPoint[2]) / 2;
				}
			}
		}

		return new FreeVertex(concaveCenter);
	}
}

float * SphereParam::FindIntersectionPoint(float *center, float *vec, float t)
{
	// Normalize direction vector
	float len = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
	vec[0] /= len;
	vec[1] /= len;
	vec[2] /= len;

	float *result = new float[3];
	result[0] = center[0] + t*vec[0];
	result[1] = center[1] + t*vec[1];
	result[2] = center[2] + t*vec[2];
	
	return result;
}

float SphereParam::FindDistance(float* coords1, float* coords2)
{
	float dist = sqrt(pow((coords2[0] - coords1[0]), 2) + pow((coords2[1] - coords1[1]), 2) + pow((coords2[2] - coords1[2]), 2));
	return dist;
}

// Move each vertex to the boundary of the sphere by simply normalizing the vector and scaling by the desired radius.
void SphereParam::GetSphereVertices(Mesh *mesh)
{
	float radius = 10;

	for(int i = 0; i < mesh->verts.size(); i++)
	{
		Vertex *v = mesh->verts[i];
		// Compute the vector
		float *sphereCoords = new float[3];
		sphereCoords[0] = v->coords[0] - centroid->coords[0];
		sphereCoords[1] = v->coords[1] - centroid->coords[1];
		sphereCoords[2] = v->coords[2] - centroid->coords[2];

		// Normalize the vector
		float len = sqrt(pow(sphereCoords[0], 2) + pow(sphereCoords[1], 2) + pow(sphereCoords[2], 2));
		sphereCoords[0] /= len;
		sphereCoords[1] /= len;
		sphereCoords[2] /= len;

		// Scale wrt radius
		sphereCoords[0] *= radius;
		sphereCoords[1] *= radius;
		sphereCoords[2] *= radius;

		mesh->verts[i]->coords = sphereCoords;
	}
}

bool SphereParam::isNegated(float *a, float *b)
{
	if(a[0] >= 0 && a[1] >= 0 && a[2] >= 0)
		if(b[0] <= 0 && b[1] <= 0 && b[2] <= 0)
			return true;
	if(a[0] >= 0 && a[1] >= 0 && a[2] <= 0)
		if(b[0] <= 0 && b[1] <= 0 && b[2] >= 0)
			return true;
	if(a[0] >= 0 && a[1] <= 0 && a[2] >= 0)
		if(b[0] <= 0 && b[1] >= 0 && b[2] <= 0)
			return true;
	if(a[0] >= 0 && a[1] <= 0 && a[2] <= 0)
		if(b[0] <= 0 && b[1] >= 0 && b[2] >= 0)
			return true;
	if(a[0] <= 0 && a[1] >= 0 && a[2] >= 0)
		if(b[0] >= 0 && b[1] <= 0 && b[2] <= 0)
			return true;
	if(a[0] <= 0 && a[1] >= 0 && a[2] <= 0)
		if(b[0] >= 0 && b[1] <= 0 && b[2] >= 0)
			return true;
	if(a[0] <= 0 && a[1] <= 0 && a[2] >= 0)
		if(b[0] >= 0 && b[1] >= 0 && b[2] <= 0)
			return true;
	if(a[0] <= 0 && a[1] <= 0 && a[2] <= 0)
		if(b[0] >= 0 && b[1] >= 0 && b[2] >= 0)
			return true;

	return false;
}

int SphereParam::FindCollisions(Mesh *mesh)
{
	//// Store the original face normals.
	//vector<float *> oldNormals;

	//for(int i = 0; i < mesh->tris.size(); i++)
	//{
	//	oldNormals.push_back(mesh->tris[i]->normal);
	//}

	//// Compute the normals after the mesh is projected onto sphere.
	//mesh->computeNormals();

	//// Store the triangles that are not correctly oriented.
	//vector<int> invertedTris;
	//for(int i = 0; i < mesh->tris.size(); i++)
	//{
	//	float *newNormal = mesh->tris[i]->normal;
	//	if(!(oldNormals[i][0] == newNormal[0] && oldNormals[i][1] == newNormal[1] && oldNormals[i][2] == newNormal[2]))
	//		int a = 0;
	//	if(isNegated(oldNormals[i], newNormal) == true)
	//	{
	//		invertedTris.push_back(i);
	//	}
	//}

	int invertedTris = 0;
	// Send rays from centroid to each vertex. If it hits another triangle, there is a collision.
	for(int i = 0; i < mesh->verts.size(); i++)
	{
		float *centerToV = new float[3];
		centerToV[0] = mesh->verts[i]->coords[0] - centroid->coords[0];
		centerToV[1] = mesh->verts[i]->coords[1] - centroid->coords[1];
		centerToV[2] = mesh->verts[i]->coords[2] - centroid->coords[2];
		int n = 0;
		for(int j = 0; j < mesh->tris.size(); j++)
		{
			float res = triangleIntersection(mesh, j, centerToV, centroid->coords);
			if(res > -1)
				n++;
		}
		if(n > 1)
			invertedTris++;
	}

	cout << "inverted triangles count: " << invertedTris << endl;
	return invertedTris;
}

Vec3f SphereParam::subtract(const Vec3f &a, const Vec3f &b)
{
	Vec3f a_minus_b;
	a_minus_b.x = a.x - b.x;
	a_minus_b.y = a.y - b.y;
	a_minus_b.z = a.z - b.z;
	return a_minus_b;
}

float SphereParam::triangleIntersection(Mesh *mesh, int tr, float* vec, float *center)
{
	Triangle *tri = mesh->tris[tr];
	Vec3f a;
	a.x = mesh->verts[tri->v1i]->coords[0];
	a.y = mesh->verts[tri->v1i]->coords[1];
	a.z = mesh->verts[tri->v1i]->coords[2];

	Vec3f b;
	b.x = mesh->verts[tri->v2i]->coords[0];
	b.y = mesh->verts[tri->v2i]->coords[1];
	b.z = mesh->verts[tri->v2i]->coords[2];

	Vec3f c;
	c.x = mesh->verts[tri->v3i]->coords[0];
	c.y = mesh->verts[tri->v3i]->coords[1];
	c.z = mesh->verts[tri->v3i]->coords[2];

	Vec3f o;
	o.x = center[0];
	o.y = center[1];
	o.z = center[2];

	float len = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
	Vec3f d;
	d.x = vec[0] / len;
	d.y = vec[1] / len;
	d.z = vec[2] / len;

	Vec3f a_minus_b = subtract(a, b);
	Vec3f a_minus_c = subtract(a, c);
	Vec3f a_minus_o = subtract(a, o);

	float detA = determinant(a_minus_b, a_minus_c, d);

	if(detA == 0.0)
	{
		return -1;
	}

	float t = (determinant(a_minus_b, a_minus_c, a_minus_o))/detA;
	if(t <= 0.0)
	{
		return -1;
	}

	float gamma = (determinant(a_minus_b,a_minus_o, d))/detA;
	if(gamma < 0 || gamma > 1)
	{
		return -1;
	}

	float beta = (determinant(a_minus_o, a_minus_c, d))/detA;
	if(beta < 0 || beta > (1 - gamma))
	{
		return -1;
	}

	return t;
}

float SphereParam::determinant(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2)
{
	return v0.x * (v1.y*v2.z - v2.y*v1.z)
			+ v0.y * (v2.x*v1.z - v1.x*v2.z)
			+ v0.z * (v1.x*v2.y - v1.y*v2.x);
}

// Move each vertex to the center of its neighbors.
void SphereParam::RelocateVertices(Mesh *mesh)
{
	for(int i = 0; i < mesh->verts.size(); i++)
	{
		Vertex *v = mesh->verts[i];
		float *center = new float[3];
		center[0] = 0;
		center[1] = 0;
		center[2] = 0;

		// Find the coordinates of the center of its neighbors.
		for(int j = 0; j < v->vertList.size(); j++)
		{
			Vertex *neighbor = mesh->verts[v->vertList[j]];
			center[0] += neighbor->coords[0];
			center[1] += neighbor->coords[1];
			center[2] += neighbor->coords[2];
		}

		center[0] /= v->vertList.size();
		center[1] /= v->vertList.size();
		center[2] /= v->vertList.size();

		mesh->verts[i]->coords = center;
	}
}

//vector<FreeVertex *> SphereParam::RelocateVertices(Mesh *mesh)
//{
//	int numberOfVertices = mesh->verts.size();
//
//	vector<FreeVertex *> resultVertices;
//	Eigen::MatrixXf W = Eigen::MatrixXf::Zero(numberOfVertices, numberOfVertices);
//	Eigen::VectorXf bx(numberOfVertices);
//	Eigen::VectorXf by(numberOfVertices);
//
//	for(int i = 0; i < numberOfVertices; i++)
//	{
//		bx(i) = 0;
//		by(i) = 0;
//
//		for(int j = 0; j < numberOfVertices; j++)
//		{
//			if(mesh->isNeighbor(i, j) == true)
//			{
//				W(i, j) = 1;
//			}
//		}
//	}
//
//	for(int i = 0; i < numberOfVertices; i++)
//	{
//		// Total value of the row
//		W(i, i) = -1 * W.rowwise().sum()(i);
//	}
//
//	// xx and xy store the coordinates of the new vertices which lie on the circle.
//	// xx and xy are ordered with respect to the global idx of the vertices.
//	Eigen::VectorXf xx = W.lu().solve(bx);
//	Eigen::VectorXf xy = W.lu().solve(by);
//
//	
//	for(int i = 0; i < numberOfVertices; i++)
//	{
//		float *coords = new float[3];
//		coords[0] = xx(i);
//		coords[1] = xy(i);
//		coords[2] = 0;
//		FreeVertex *v = new FreeVertex(coords);
//		resultVertices.push_back(v);
//	}
//
//	return resultVertices;
//}


void SphereParam::SphericalParameterization(Mesh *mesh)
{
	// Find centroid of the mesh.
	centroid = FindCentroid(mesh);

	cout << "centroid" << endl;
	cout << centroid->coords[0] << endl;
	cout << centroid->coords[1] << endl;
	cout << centroid->coords[2] << endl;

	GetSphereVertices(mesh);
	/*int n = FindCollisions(mesh);

	while(FindCollisions(mesh) > 0)
	{
		RelocateVertices(mesh);
		GetSphereVertices(mesh);
	}*/
	// Find sphere points by sending rays from centroid to every vertex.
	/*for(int i = 0; i < 45; i++)
	{
		RelocateVertices(mesh);
		GetSphereVertices(mesh);
	}*/

	//FindCollisions(mesh);
}
