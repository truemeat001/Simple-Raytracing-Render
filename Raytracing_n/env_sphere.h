#pragma once
#ifndef _HITABLE_ENVSPHERE_H_
#define _HITABLE_ENVSPHERE_H_
#include "onb.h"
#include "hitable.h"

class env_sphere :public hitable
{
public:
	env_sphere();
	env_sphere(vec3 cen, float r, material *m):center(cen), radius(r), mat_ptr(m){}
	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec, bool is_medium = false) const;
	virtual bool bounding_box(float t0, float t1, aabb& box)const;
	virtual float pdf_value(const vec3& o, const vec3& v) const;

	virtual vec3 random(const vec3& o) const;
	vec3 center;
	float radius;
	material* mat_ptr;
};

bool env_sphere::bounding_box(float t0, float t1, aabb& box)const {
	box = aabb(center - vec3(radius, radius, radius), center + vec3(radius, radius, radius));
	return true;
}

bool env_sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec, bool is_medium) const {
	vec3 oc = center - r.origin();
	float proj = dot(r.direction(), oc);
	float distance2 = oc.squared_length() - proj * proj;
	float a = sqrt(radius * radius - distance2);
	rec.t = a + proj;
	rec.p = r.point_at_parameter(rec.t);
	get_sphere_uv((rec.p - center) / radius, rec.u, rec.v);
	rec.normal = (center - rec.p) / radius;
	rec.mat_ptr = mat_ptr;
	return true;
}

vec3 env_sphere::random(const vec3& o) const {
	vec3 direction = o - center;
	float distance_squred = direction.squared_length();
	onb uvw;
	uvw.build_from_w(direction);
	return uvw.local(random_to_sphere(radius, distance_squred));
}

float env_sphere::pdf_value(const vec3& o, const vec3& v) const {
	return 1;
}


#endif	// _HITABLE_ENVSPHERE_H_