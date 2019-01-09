
#ifndef MODEL_H
#define MODEL_H
#include "common.h"
#include "geometry.h"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>

class geometry;
class material;

class model
{
public:
	model(const string& filename, bool flipUVs, bool flipWindingOrder, material* mat, vec3 scale);
	~model();


	const std::vector<geometry*>& geometries() const;
	hitable** genhitablemodel();
	int gettrianglecount();
private:
	std::vector<geometry*> geometries_;
};

model::model(const string& filename, bool flipUVs, bool flipWindingOrder, material *mat, vec3 scale)
{
	Assimp::Importer impoter;

	int flags = aiProcess_Triangulate | aiProcess_JoinIdenticalVertices | aiProcess_SortByPType;// ;

	if (flipUVs)
	{
		flags |= aiProcess_FlipUVs;
	}
	if (flipWindingOrder)
	{
		flags |= aiProcess_FlipWindingOrder;
	}
	
	const aiScene* scene = impoter.ReadFile(filename, flags);
	if (scene == nullptr)
	{

	}

	if (scene->HasMeshes())
	{
		for (int i = 0; i < scene->mNumMeshes; i++)
		{
			geometry *geo = new geometry(*(scene->mMeshes[i]), mat, scale);
			geometries_.push_back(geo);
		}
	}
	
}


model::~model()
{
	for (geometry* geo : geometries_)
	{
		delete geo;
	}
}

const std::vector<geometry*>& model::geometries() const
{
	return geometries_;
}

hitable** model::genhitablemodel() 
{
	return geometries_[0]->gethitablegeometry();
}

int model::gettrianglecount()
{
	return geometries_[0]->gettrianglecount();
}
#endif

