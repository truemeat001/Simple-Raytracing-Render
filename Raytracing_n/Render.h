#pragma once
#ifndef __CRENDER_H__
#define __CRENDER_H__
#include "common.h"
#include <ctime>
#include <chrono>
#include <thread>
#include <mutex>
#include "stdlib.h"

class ray;
class hitable;
class camera;
class vec3;

class CRender
{
public:
	CRender();
	CRender(int maxDepth, int width, int height, int num_ray, int thread_count, char* filePath);
	~CRender();

	void Run();
private:
	std::ofstream outfile;
	int **colors;
	int curfinishedcount;
	std::mutex g_lock;
	std::vector<std::thread> threads;
	char* filePath;


	int maxDepth;
	int nx;
	int ny;
	int ns;
	int thread_count;
	int finishedPixel;
	bool isfinished;

	int getpixel();
	vec3 color(const ray& r, hitable *world, hitable *light_shape, int *depth);
	void cornell_box(hitable **scene, camera **cam, hitable **hlist, float aspect);
	void soldier_scene(hitable **scene, camera **cam, hitable **hlist, float aspect);
	void renderthread(hitable* world, hitable* hlist, camera *cam, double** samplepoints);
	double **sobol_points(unsigned N, unsigned D, const char *dir_file);

	vec3 de_nan(const vec3& c);
	
	
};



#endif