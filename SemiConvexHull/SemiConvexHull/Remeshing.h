#pragma once

#include <Eigen/Dense>
#include "Mesh.h"
#include "Sampling.h"
#include <map>
#include <set>
  
#define PI 3.1415926535 

class Remeshing
{
private:
	float FindDistance(float *a, float *b);
	float* FindCrossProduct(float *a, float *b);
	float FindDotProduct(float *a, float *b);
	float FindLength(float *a);
	int GetIndex(const vector<int> list, int v);
	float FindAvgAdjLength(Mesh *mesh, int v);
	float* FindTriangleCenter(Mesh *mesh, int t);
	
	void FlipEdge(Mesh *mesh, int v1, int v2, float *normal);
	float FindMinAngleOfTriangle(float *v1coords, float *v2coords, float *v3coords);
	bool IsValidEdgeFlip(Mesh *mesh, int v1, int v2);
	void RemoveTriangleInfoFromVertices(Mesh *mesh, int id);
	void AddTriangleInfoToVertices(Mesh *mesh, int id);
	void RemoveEdgeInfoFromVertices(Mesh *mesh, int edgeId);
	void AddEdgeInfoToVertices(Mesh *mesh, int edgeId);
	bool IsConvexQuad(Mesh *mesh, int otherVertex, int v1, int v2, int v3);
	float* FindPlaneNormal(float *a, float *b, float *c);
public:
	int FindNN(Mesh *mesh, int v, vector<float *> samples);
	void Remesh(Mesh *mesh, float avgLen);
	void OptimizeMesh(Mesh *mesh, float avgLen, vector<float *> samples);
	void ShiftVertices(Mesh *mesh, float gamma, vector<vector<float>> dE);
	void ScaleVertices(Mesh *mesh, float scale);
	vector<vector<float>> CalculateEnergyDerivative(Mesh *mesh, float w, float avgEdgeLen, vector<float *> samples);
};
