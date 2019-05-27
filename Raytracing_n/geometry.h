#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "common.h"
#include "triangle.h"
#include "material.h"
#include "bvh.h"
#include <assimp/scene.h>

class geometry {
public :
	geometry(aiMesh& mesh, material *mat, vec3 scale);
	hitable** gethitablegeometry();
	int gettrianglecount();
private:
	std::vector<vec3> vertices_;
	std::vector<vec3> texcoords_;
	std::vector<triangle*> triangles;
	int triangleCount;
	material* mp;
};

geometry::geometry(aiMesh& mesh, material *mat, vec3 scale)
{
	mp = mat;
	vertices_.reserve(mesh.mNumVertices);
	for (int i = 0; i < mesh.mNumVertices; i++)
	{
		vertices_.push_back(vec3(mesh.mVertices[i].x, mesh.mVertices[i].y, mesh.mVertices[i].z));
	}
	
	int uvChannelCount = mesh.GetNumUVChannels();
	if (uvChannelCount > 0)
	{
		texcoords_.reserve(mesh.mNumVertices);
		for (int i = 0; i < uvChannelCount; i++)
		{
			aiVector3D* aiTextureCoordinates = mesh.mTextureCoords[i];
			for (int j = 0; j < mesh.mNumVertices; j++)
			{
				texcoords_.push_back(vec3(reinterpret_cast<const float*>(&aiTextureCoordinates[j])));
			}
		}
	}
	

	if (mesh.HasFaces())
	{
		triangleCount = mesh.mNumFaces;
		triangles.reserve(triangleCount);
		aiVector3D *normals = mesh.mNormals;
		vec3 p[3];
		vec3 uv[3];
		for (int i = 0; i < triangleCount; i++)
		{
			aiFace* face = &mesh.mFaces[i];
			for (int j = 0; j < 3; j++)
			{
				p[j] = vertices_[face->mIndices[j]];
				p[j] = vec3(p[j].x() * scale.x(), p[j].y() * scale.y(), p[j].z() * scale.z()); 
				if (uvChannelCount > 0)
				{
					uv[j] = texcoords_[face->mIndices[j]];
				}
			}
			
			triangle* tri = new triangle(p[0], p[1], p[2],  mp, uv[0], uv[1], uv[2]);
			triangles.push_back(tri);
		}
	}
}

hitable** geometry::gethitablegeometry()
{
	hitable** trianglelist = new hitable*[triangleCount];
	for (int i = 0; i < triangleCount; i++)
	{
		trianglelist[i] = triangles[i];
	}
	return trianglelist;// new bvh_node(trianglelist, triangleCount, 0, 1);
}

int geometry::gettrianglecount()
{
	return triangleCount;
}

inline vec3 SphericalDirection(float sinTheta, float cosTheta, float phi) {
	return vec3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
}
//
//template <typename T>
//class Point2 {
//public:
//	explicit Point2(const Point3<T> &p)
//};
//
//template<typename T>
//class Point3 {
//public:
//	Point3() { x = y = z = 0; }
//	Point3(T x, T y, T z): x(x), y(y), z(z){}
//	template <typename U>
//	explicit Point3(const Point3<U> &p)
//		:x((T)p.x, y((T)p.y), z((T)p.z)) {
//	}
//
//	Point3<T> &operator+=(const Point3<T> &p) {
//		x += p.x;
//		y += p.y;
//		z += p.z;
//	}
//
//	T x, y, z;
//};
//
//typedef Point2<float> Point2f;

#endif