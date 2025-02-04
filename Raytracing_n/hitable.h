#pragma once
#ifndef HITABLEH
#define HITABLEH

#include "ray.h"
#include "aabb.h"
#include "mathf.h"
class material;

void get_sphere_uv(const vec3& p, float& u, float& v) {
	float phi = atan2(p.z(), p.x());
	float theta = asin(p.y());
	u = 1 - (phi + M_PI) / (2 * M_PI);
	v = (theta + M_PI / 2) / M_PI;
}

struct hit_record
{
	float t;
	float u;
	float v;
	vec3 p;
	vec3 normal;
	material *mat_ptr;
};

class hitable {
public:
	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec, bool is_medium = false) const = 0;
	virtual bool bounding_box(float t0, float t1, aabb& box) const = 0;
	virtual float pdf_value(const vec3& o, const vec3& direction) const { return 0.0; };
	virtual vec3 random(const vec3& o) const { return vec3(1, 0, 0); };
};

class translate : public hitable {
public :
	translate(hitable *p, const vec3& displacement) :ptr(p), offset(displacement) {}
	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec, bool is_medium = false) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const;
	hitable *ptr;
	vec3 offset;
};

bool translate::hit(const ray& r, float t_min, float t_max, hit_record& rec, bool is_medium) const {
	ray moved_r(r.origin() - offset, r.direction(), r.time());
	if (ptr->hit(moved_r, t_min, t_max, rec, is_medium)) {
		rec.p += offset;
		return true;
	}
	else
		return false;
}

bool translate::bounding_box(float t0, float t1, aabb& box) const {
	if (ptr->bounding_box(t0, t1, box)) {
		box = aabb(box.Max() + offset, box.Max() + offset);
		return true;
	}
	else
		return false;
}



class rotate_y : public hitable {
public:
	rotate_y(hitable* p, float angle);
	virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec, bool is_medium = false) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const {
		box = bbox; 
		return hasbox;
	}

	hitable *ptr;
	float sin_theta;
	float cos_theta;
	bool hasbox;
	aabb bbox;
};

rotate_y::rotate_y(hitable *p, float angle) :ptr(p) {
	float radians = (M_PI / 180.) * angle;
	sin_theta = sin(radians);
	cos_theta = cos(radians);
	hasbox = ptr->bounding_box(0, 1, bbox);
	vec3 _min(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 _max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	for(int i =0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
			{
				float x = i * bbox.Max().x() + (1 - i)*bbox.Min().x();
				float y = j * bbox.Max().y() + (1 - j)*bbox.Min().y();
				float z = k * bbox.Max().z() + (1 - k)*bbox.Min().z();
				float newx = cos_theta * x + sin_theta * z;
				float newz = -sin_theta * x + cos_theta * z;
				vec3 tester(newx, y, newz);
				for (int c = 0; c < 3; c++)
				{
					if (tester[c] > _max[c])
						_max[c] = tester[c];
					if (tester[c] < _min[c])
						_min[c] = tester[c];
				}
			}
	bbox = aabb(_min, _max);
}

bool rotate_y::hit(const ray& r, float t_min, float t_max, hit_record& rec, bool is_medium ) const {
	vec3 origin = r.origin();
	vec3 direction = r.direction();
	origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[2];
	origin[2] = sin_theta * r.origin()[0] + cos_theta * r.origin()[2];
	direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[2];
	direction[2] = sin_theta * r.direction()[0] + cos_theta * r.direction()[2];
	ray rotated_r(origin, direction, r.time());
	if (ptr->hit(rotated_r, t_min, t_max, rec, is_medium))
	{
		vec3 p = rec.p;
		vec3 normal = rec.normal;
		p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
		p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];
		normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
		normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];
		rec.p = p;
		rec.normal = normal;
		return true;
	}
	else
		return false;

}


class rotate_x : public hitable {
public:
	rotate_x(hitable* p, float angle);
	virtual bool hit(const ray &r, float t_min, float t_max, hit_record& rec, bool is_medium ) const;
	virtual bool bounding_box(float t0, float t1, aabb& box) const {
		box = bbox;
		return hasbox;
	}

	hitable* ptr;
	float sin_theta;
	float cos_theta;
	bool hasbox;
	aabb bbox;
};

rotate_x::rotate_x(hitable *p, float angle): ptr(p) {
	float radians = (M_PI / 180.) * angle;
	sin_theta = sin(radians);
	cos_theta = cos(radians);
	hasbox = ptr->bounding_box(0, 1, bbox);
	vec3 _min(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 _max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 2; k++)
			{
				float x = i * bbox.Max().x() + (1 - i)*bbox.Min().x();
				float y = j * bbox.Max().y() + (1 - j)*bbox.Min().y();
				float z = k * bbox.Max().z() + (1 - k)*bbox.Min().z();
				float newy = cos_theta * y + sin_theta * z;
				float newz = -sin_theta * y + cos_theta * z;
				vec3 tester(x, newy, newz);
				for (int c = 0; c < 3; c++)
				{
					if (tester[c] > _max[c])
						_max[c] = tester[c];
					if (tester[c] < _min[c])
						_min[c] = tester[c];
				}

			}
	bbox = aabb(_min, _max);
}

bool rotate_x::hit(const ray &r, float t_min, float t_max, hit_record& rec, bool is_medium = false) const
{
	vec3 origin = r.origin();
	vec3 direction = r.direction();
	origin[1] = cos_theta * r.origin()[1] - sin_theta * r.origin()[2];
	origin[2] = sin_theta * r.origin()[1] + cos_theta * r.origin()[2];
	direction[1] = cos_theta * r.direction()[1] - sin_theta * r.direction()[2];
	direction[2] = sin_theta * r.direction()[1] + cos_theta * r.direction()[2];
	ray rotated_r(origin, direction, r.time());
	if (ptr->hit(rotated_r, t_min, t_max, rec, is_medium))
	{
		vec3 p = rec.p;
		vec3 normal = rec.normal;
		p[1] = cos_theta * rec.p[1] + sin_theta * rec.p[2];
		p[2] = -sin_theta * rec.p[1] + cos_theta * rec.p[2];
		normal[1] = cos_theta * rec.normal[1] + sin_theta * rec.normal[2];
		normal[2] = -sin_theta * rec.normal[1] + cos_theta * rec.normal[2];
		rec.p = p;
		rec.normal = normal;
		return true;
	}
	else
		return false;
}

#endif




