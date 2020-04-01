#include "SphereParam.h"

FreeVertex* SphereParam::FindCentroid(Mesh *mesh)
{
	// Find the average of all vertex coordinates.
	float *center = new float[3];
	center[0] = 0;
	center[1] = 0;
	center[2] = 0;

	// burayi commentle asagiyi dene
	for(int i = 0; i < mesh->verts.size(); i++)
	{
		center[0] += mesh->verts[i]->coords[0];
		center[1] += mesh->verts[i]->coords[1];
		center[2] += mesh->verts[i]->coords[2];
	}

	center[0] /= mesh->verts.size();
	center[1] /= mesh->verts.size();
	center[2] /= mesh->verts.size();
	
	// Find an outside point by averaging farthest points.
	// Sampling fps = new Sampling(5);
	// vector<int> samples = fps->FPS(mesh);
	//
	// for(int i = 0; i < samples.size(); i++)
	// {
	//	center[0] += mesh->verts[samples[i]]->coords[0];
	//	center[1] += mesh->verts[samples[i]]->coords[1];
	//	center[2] += mesh->verts[samples[i]]->coords[2];
	// }
	
	// center[0] /= samples.size();
	// center[1] /= samples.size();
	// center[2] /= samples.size();

	if(isConvex == true)
		return new FreeVertex(center);
	else
	{
		// Find the closest vertex to found center.
		float minDist = numeric_limits<float>::infinity();;
		int minV;
		for(int i = 0; i < mesh->verts.size(); i++)
		{
			float dist = FindDistance(center, mesh->verts[i]->coords);
			if(dist < minDist)
			{
				minDist = dist;
				minV = i;
			}
		}

		float *closestCoords = mesh->verts[minV]->coords;
		float *ray = new float[3];
		ray[0] = closestCoords[0] - center[0];
		ray[1] = closestCoords[1] - center[1];
		ray[2] = closestCoords[2] - center[2];

		float *concaveCenter = new float[3];

		// Triangle intersection
		for(int i = 0; i < mesh->tris.size(); i++)
		{
			if(mesh->tris[i]->v1i != minV && mesh->tris[i]->v2i != minV && mesh->tris[i]->v3i != minV)
			{
				float hitResult = TriangleIntersection(mesh, i, ray, center);
				if(hitResult > 0)
				{
					// Find intersection point on the triangle.
					float *hitPoint = FindIntersectionPoint(center, ray, hitResult);
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
	// Normalize direction vector(vec)
	float len = FindVectorLen(vec);
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

float SphereParam::FindVectorLen(float *vec)
{
	return sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
}

// Move each vertex to the boundary of the sphere by simply normalizing the vector and scaling by the desired radius.
void SphereParam::GetSphereVertices(Mesh *mesh)
{
	float radius = 10;

	// Find sphere points by sending rays from centroid to every vertex.
	for(int i = 0; i < mesh->verts.size(); i++)
	{
		Vertex *v = mesh->verts[i];
		// Compute the vector
		float *sphereCoords = new float[3];
		sphereCoords[0] = v->coords[0] - centroid->coords[0];
		sphereCoords[1] = v->coords[1] - centroid->coords[1];
		sphereCoords[2] = v->coords[2] - centroid->coords[2];

		// Normalize the vector
		float len = FindVectorLen(sphereCoords);
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
	
	// Find the centers of triangles after mapping to sphere
	vector<float *> centers;
	for(int i = 0; i < mesh->tris.size(); i++)
	{
		float *v1Coords = mesh->tris[i]->v1i->coords;
		float *v2Coords = mesh->tris[i]->v2i->coords;
		float *v3Coords = mesh->tris[i]->v3i->coords;
		
		float *center = new float[3];
		center[0] = (v1Coords[0] + v2Coords[0] + v3Coords[0]) / 3;
		center[1] = (v1Coords[1] + v2Coords[1] + v3Coords[1]) / 3;
		center[2] = (v1Coords[2] + v2Coords[2] + v3Coords[2]) / 3;
		centers.push_back(center);
	}
	
	int invertedTris = 0;
	for(int i = 0; i < mesh->tris.size(); i++)
	{
		if(centers[i][0] * mesh->tris[i]->normal[0] 
		  + centers[i][1] * mesh->tris[i]->normal[1] 
		  + centers[i][2] * mesh->tris[i]->normal[2] < 0)
		{
			invertedTris++;
		}
	}

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
		// Keep neighboring triangles of the vertex.
		vector<int> neighborTris = mesh->verts[i]->triList;
		
		float *centerToV = new float[3];
		centerToV[0] = mesh->verts[i]->coords[0] - centroid->coords[0];
		centerToV[1] = mesh->verts[i]->coords[1] - centroid->coords[1];
		centerToV[2] = mesh->verts[i]->coords[2] - centroid->coords[2];
		
		for(int j = 0; j < mesh->tris.size(); j++)
		{
			// If the selected vertex is not part of the triangle, try intersection with the triangle.
			if(find(neigborTris.begin(), neighborTris.end(), j) != neighborTris.end())
			{
				float hitResult = TriangleIntersection(mesh, j, centerToV, centroid->coords);
				if(hitResult > 0)
					invertedTris++;
			}
		}
	}

	cout << "inverted triangles count: " << invertedTris << endl;
	return invertedTris;
}

float* SphereParam::Subtract(float *a, float *b)
{
	float *result = new float[3];
	result[0] = a[0] - b[0];
	result[1] = a[1] - b[1];
	result[2] = a[2] - b[2];
	return result;
}

float SphereParam::TriangleIntersection(Mesh *mesh, int tr, float* vec, float *center)
{
	Triangle *tri = mesh->tris[tr];
	float *v1Coords = mesh->verts[tri->v1i]->coords;
	float *v2Coords = mesh->verts[tri->v2i]->coords;
	float *v3Coords = mesh->verts[tri->v3i]->coords;

	// Normalize direction vector(vec)
	float len = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
	vec[0] /= len;
	vec[1] /= len;
	vec[2] /= len;

	float *a_minus_b = Subtract(v1Coords, v2Coords);
	float *a_minus_c = Subtract(v1Coords, v3Coords);
	float *a_minus_o = Subtract(v1Coords, center);

	float detA = Determinant(a_minus_b, a_minus_c, vec);

	if(detA == 0.0)
	{
		return -1;
	}

	float t = (Determinant(a_minus_b, a_minus_c, a_minus_o))/detA;
	if(t <= 0.0)
	{
		return -1;
	}

	float gamma = (Determinant(a_minus_b,a_minus_o, vec))/detA;
	if(gamma < 0 || gamma > 1)
	{
		return -1;
	}

	float beta = (Determinant(a_minus_o, a_minus_c, vec))/detA;
	if(beta < 0 || beta > (1 - gamma))
	{
		return -1;
	}

	return t;
}

float SphereParam::Determinant(float *v0, float *v1, float *v2)
{
	return v0[0] * (v1[1] * v2[2] - v2[1] * v1[2])
		+ v0[1] * (v2[0] * v1[2] - v1[0] * v2[2])
		+ v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]);
}

// Move each vertex to the center of its neighbors.
void SphereParam::RelocateVertices(Mesh *mesh)
{
	for(int i = 0; i < mesh->verts.size(); i++)
	{
		vector<int> neighborVertices = mesh->verts[i]->vertList;
		float *center = new float[3];
		center[0] = 0;
		center[1] = 0;
		center[2] = 0;

		// Find the coordinates of the center of its neighbors.
		for(int j = 0; j < neighborVertices.size(); j++)
		{
			Vertex *neighbor = mesh->verts[neighborVertices[j]];
			center[0] += neighbor->coords[0];
			center[1] += neighbor->coords[1];
			center[2] += neighbor->coords[2];
		}

		center[0] /= neighborVertices.size();
		center[1] /= neighborVertices.size();
		center[2] /= neighborVertices.size();

		mesh->verts[i]->coords = center;
	}
}

//void SphereParam::RelocateVertices(Mesh *mesh)
//{
//	int numberOfVertices = mesh->verts.size();
//
//	Eigen::MatrixXf W = Eigen::MatrixXf::Zero(numberOfVertices, numberOfVertices);
//	Eigen::VectorXf bx(numberOfVertices);
//	Eigen::VectorXf by(numberOfVertices);
//	Eigen::VectorXf bz(numberOfVertices);
//
//	for(int i = 0; i < numberOfVertices; i++)
//	{
//		if(i == 0)
//		{
//			bx(i) = mesh->verts[i]->coords[0];
//			by(i) = mesh->verts[i]->coords[1];
//			bz(i) = mesh->verts[i]->coords[2];
//		}
//		else
//		{
//			bx(i) = 0;
//			by(i) = 0;
//			bz(i) = 0;
//		}
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
//	Eigen::VectorXf xz = W.lu().solve(bz);
//
//	
//	for(int i = 0; i < numberOfVertices; i++)
//	{
//		float *coords = new float[3];
//		coords[0] = xx(i);
//		coords[1] = xy(i);
//		coords[2] = xz(i);
//		mesh->verts[i]->coords = coords;
//	}
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
	
	/*int invertedTris = FindCollisions(mesh);
	while(invertedTris > 0)
	{
		cout << "Number of inverted triangles: " << invertedTris << endl;
		cout << "Relocating vertices..." << endl;
		RelocateVertices(mesh);
		GetSphereVertices(mesh);
		// Recalculate the number of inverted triangles.
		invertedTris = FindCollisions(mesh);
	}*/
	
	/*for(int i = 0; i < 45; i++)
	{
		RelocateVertices(mesh);
		GetSphereVertices(mesh);
	}*/

	//FindCollisions(mesh);
}
