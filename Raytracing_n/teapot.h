#pragma once
#ifndef _HITABLE_TEAPOT_H_
#define _HITABLE_TEAPOT_H_
#include "common.h"
#include "geometry.h"
#include "teapotdata.h"

class geometry;

class teapot
{
public:
	teapot(float scale, material *mat) : scale(scale), mat(mat) {}
	float scale;
	std::vector<triangle*> triangles;
	int triangleCount;
	material *mat;
	// Compute the position of a point along a Bezier Curve at t [0:1]
	vec3 evalBezierCurve(const vec3 *p, const float &t)
	{
		float b0 = (1 - t)*(1 - t)*(1 - t);
		float b1 = 3 * t*(1 - t)*(1 - t);
		float b2 = 3 * t*t*(1 - t);
		float b3 = t * t*t;
		
		return p[0] * b0 + p[1] * b1 + p[2] * b2 + p[3] * b3;
	}

	vec3 evalBezierPatch(const vec3* controlPoints, const float &u, const float &v)
	{
		vec3 uCurve[4];
		for (size_t i = 0; i < 4; i++)
		{
			uCurve[i] = evalBezierCurve(controlPoints + 4 * i, u);
		}
		return evalBezierCurve(uCurve, v);
	}

	vec3 derivBezier(const vec3 *p, const float &t)
	{
		return -3 * (1 - t)*(1 - t) * p[0] +
			(3 * (1 - t)*(1 - t) - 6 * t*(1 - t))*p[1] +
			(6 * t* (1 - t) - 3 * t*t) * p[2] +
			3 * t*t*p[3];
	}

	// Compute the derivative of a point on Bezier patch along the u parametric direction
	vec3 dUBezier(const vec3 * controlPoints, const float &u, const float &v)
	{
		vec3 p[4];
		vec3 vCurve[4];
		for(int i = 0; i < 4; i++)
		{
			p[0] = controlPoints[i];
			p[1] = controlPoints[4 + i];
			p[2] = controlPoints[8 + i];
			p[3] = controlPoints[12 + i];
		}

		return derivBezier(vCurve, u);
	}

	// Compute the derivative of a point on Bezier patch along the v parametric direction
	vec3 dVBezier(const vec3 *controlPoints, const float &u, const float &v)
	{
		vec3 uCurve[4];
		for (int i = 0; i < 4; i++)
		{
			uCurve[i] = evalBezierCurve(controlPoints + 4 * i, u);
		}

		return derivBezier(uCurve, v);
	}

	// Generate a poly-mesh Utah teapot out of a Bezier patches
	hitable** createPloyTeapot() {
		int divs = 100;
		vec3* P = new vec3[(divs + 1)* (divs + 1)];
		int* nvertices = new int[divs * divs];
		int* vertices = new int[divs * divs * 4];
		vec3* n = new vec3[(divs + 1)* (divs + 1)];
		vec3* st = new vec3[(divs + 1) * (divs + 1)];

		// face connectivity - all patches are subdivided the same way so there
		// share the same topplogy and uvs
		for (int j = 0, k = 0; j < divs; ++j)
		{
			for (int i = 0; i < divs; ++i, ++k) {
				nvertices[k] = 4;
				vertices[k * 4] = (divs + 1) * j + i;
				vertices[k * 4 + 1] = (divs + 1) * j + i + 1;
				vertices[k * 4 + 2] = (divs + 1) * (j + 1) + i + 1;
				vertices[k * 4 + 3] = (divs + 1) * (j + 1) + i;
			}
		}

		vec3 controlPoints[16];
		for (int np = 0; np < kTeapotNumPatches; ++np) {		// kTeapotNumPatches
			for (int i = 0; i < 16; ++i)
			{
				controlPoints[i][0] = teapotVertices[teapotPatches[np][i] - 1][0] * scale;
				controlPoints[i][1] = teapotVertices[teapotPatches[np][i] - 1][1] * scale;
				controlPoints[i][2] = teapotVertices[teapotPatches[np][i] - 1][2] * scale;
			}
			// gereate grid
			for (int j = 0, k = 0; j <= divs; ++j)
			{
				float v = (float)j / (float)divs;
				for (int i = 0; i <= divs; ++i, ++k)
				{
					float u = (float)i / (float)divs;
					P[k] = evalBezierPatch(controlPoints, u, v);
					/*vec3 dU = dUBezier(controlPoints, u, v);
					vec3 dV = dVBezier(controlPoints, u, v);
					st[k].x = u;
					st[k].y = v;*/
				}
			}
			int numTris = 0;
			for (int i = 0; i < divs * divs; i++)
			{
				numTris += nvertices[i] - 2;

			}
			for (int i = 0, k = 0, l = 0; i < divs * divs; ++i)
			{
				for (int j = 0; j < nvertices[i] - 2; ++j)
				{
					triangle* tri = new triangle(P[vertices[k]], P[vertices[k + j + 1]], P[vertices[k + j + 2]], mat);
					triangles.push_back(tri);
					l += 3;
				}
				k += nvertices[i];
			}
			triangleCount += numTris;

		}
		
		

		//triangleCount = divs * divs * 2;
	/*	for (int j = 0; j < divs; ++j) {
			for (int i = 0; i < divs; ++i)
			{
				triangle* tri0 = new triangle(P[j*divs + i], P[j*divs + i + 1], P[(j + 1)*divs + i + 1], mat);
				triangle* tri1 = new triangle(P[j*divs + i], P[(j + 1)*divs + i + 1], P[(j + 1)*divs + i], mat);
				triangles.push_back(tri0);
				triangles.push_back(tri1);
			}
		}*/
		//for (int i = 0; i < divs * divs; i++)
		//{
		//	triangle* tri0 = new triangle(P[vertices[4 * i]], P[vertices[4 * i + 2]], P[vertices[4 * i] + 1], mat);
		//	//triangle* tri1 = new triangle(P[vertices[4 * i]], P[vertices[4 * i + 2]], P[vertices[4 * i] + 1], mat);
		//	triangle* tri1 = new triangle(P[vertices[4 * i + 2]], P[vertices[4 * i + 3]], P[vertices[4 * i] + 0], mat);
		//	triangles.push_back(tri0);
		//	triangles.push_back(tri1);
		//}

		hitable** list = new hitable*[triangleCount];
		for (int i = 0; i < triangleCount; i++)
		{
			list[i] = triangles[i];
		}
		return list;
	}

	int getTriangleCount()
	{
		return triangleCount;
	}
};

#endif  // _HITABLE_TEAPOT_H_