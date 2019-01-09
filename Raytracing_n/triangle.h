#ifndef GEOMETRY_TRIANLE_H
#define GEOMETRY_TRIANLE_H

#include "hitable.h"
#include "vec3.h"
#include <iostream>

class triangle : public hitable
{
public :
	triangle() {}
	triangle(vec3 _p0, vec3 _p1, vec3 _p2, material *mat) :
		p0(_p0), p1(_p1), p2(_p2), mp(mat) {
		points[0] = p0;
		points[1] = p1;
		points[2] = p2;
		normal = cross(p1 - p0, p2 - p0);
		normal = unit_vector(normal);
	};
	virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
	virtual bool bounding_box(float t0, float t1, aabb& box)const {
		float xmin = ffmin(p0.x(), p1.x());
		xmin = ffmin(xmin, p2.x());
		float xmax = ffmax(p0.x(), p1.x());
		xmax = ffmax(xmax, p2.x());
		float ymin = ffmin(p0.y(), p1.y());
		ymin = ffmin(ymin, p2.y());
		float ymax = ffmax(p0.y(), p1.y());
		ymax = ffmax(ymax, p2.y());
		float zmin = ffmin(p0.z(), p1.z());
		zmin = ffmin(zmin, p2.z());
		float zmax = ffmax(p0.z(), p1.z());
		zmax = ffmax(zmax, p2.z());
		box = aabb(vec3(xmin, ymin, zmin), vec3(xmax, ymax, zmax));
		return true;
	}

	virtual float pdf_value(const vec3& o, const vec3& v) const {
		hit_record rec;
		if (this->hit(ray(o, v), 0.001, FLT_MAX, rec)) {
			vec3 v01 = p1 - p0;
			vec3 v01n = v01 / v01.length();
			vec3 v02 = p2 - p0;
			vec3 v02n = v02 / v02.length();
			float cos102 = dot(v01n, v02n);
			float sin102 = sqrt(1 - cos102 * cos102);
			float h = v02.length() * sin102;
			float area = 0.5 * v01.length() * h;
			float distance_square = rec.t*rec.t*v.squared_length();
			float cosine = fabs(dot(v, rec.normal)) / v.length();
			return distance_square / (cosine * area);
		}
		else
			return 0;
	}

	virtual vec3 random(const vec3& o) const {
		float u = drand48();
		float v = drand48() * (1 - u);
		vec3 random_point = p0 * (1 - u - v) + p1 * u + p2 * v;
		return random_point - o;
	}
	// clock-wise
	material *mp;
	vec3 normal;
	vec3 p0, p1, p2;

	vec3 points[3];
};
/*
bool triangle::hit(const ray& r, float t0, float t1, hit_record& rec) const {
	vec3 v0v1 = p1 - p0;
	vec3 v0v2 = p2 - p0;

	vec3 N = cross(v0v1, v0v2);

	float area = N.length();

	float NdotRayDirection = dot(N, r.direction());
	if (fabs(NdotRayDirection) < 0.001)
		return false;

	float d = dot(N, p0);

	float t = (dot(N, r.origin()) + d) / NdotRayDirection;

	if (t < 0)
		return false;
	vec3 p = r.origin() + t * r.direction();
	vec3 c;
	vec3 edge0 = p1 - p0;
	vec3 vp0 = p - p0;
	c = cross(edge0, vp0);
	if (dot(N, c) < 0)
		return false;

	vec3 edge1 = p2 - p1;
	vec3 vp1 = p - p1;
	c = cross(edge1, vp1);
	if (dot(N, c) < 0)
		return false;

	vec3 edge2 = p0 - p2;
	vec3 vp2 = p - p2;
	c = cross(edge2, vp2);
	if (dot(N, c) < 0)
		return false;

	rec.mat_ptr = mp;
	rec.normal = normal;
	rec.p = p;
	rec.t = (p - r.origin()).length();
	return true;
}*/

bool triangle::hit(const ray& r, float t0, float t1, hit_record& rec) const {
	float u, v, t;

	vec3 e1 = p1 - p0;
	vec3 e2 = p2 - p0;

	vec3 dir = r.direction() / r.direction().length();

	vec3 P = cross(dir, e2);

	float det = dot(e1, P);
	

	vec3 T;
	if (det > 0)
	{
		T = r.origin() - p0;
	}
	else
	{
		T = p0 - r.origin();
		det = -det;
	}

	// parallel
	if (det < 0.0001)
		return false;

	// calculate u and make sure u <=1
	u = dot(T, P);
	if (u < 0.0f || u > det)
		return false;

	// Q
	vec3 Q = cross(T, e1);

	// calculate v and make u + v <= 1
	v = dot(dir, Q);
	if (v < 0.0f || v + u > det)
		return false;

	t = dot(e2, Q);

	float fInvDet = 1.0f / det;
	t *= fInvDet;
	u *= fInvDet;
	v *= fInvDet;
	if (t < 0.0001)
		return false;

	rec.u = 0;
	rec.v = 0;
	rec.mat_ptr = mp;
	rec.normal = normal;
	rec.p = (1 - u - v) * p0 + u * p1 + v * p2;
	t = (rec.p - r.origin()).length();
	rec.t = t;
	return true;
}

#endif
