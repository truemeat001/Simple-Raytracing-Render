#pragma once
#ifndef  AABBH
#define AABBH

#include "vec3.h"
#include "ray.h"
inline float ffmin(float a, float b) { return a < b ? a : b; }
inline float ffmax(float a, float b) { return a > b ? a : b; }
using namespace std;
class aabb {
public:
	aabb() {}
	aabb(const vec3& a, const vec3& b) { _min = a; _max = b; }

	vec3 Min() const { return _min; }
	vec3 Max() const { return _max; }

	/*bool hit(const ray& r, float tmin, float tmax) const {
		for (int a = 0; a < 3; a++)
		{
			float t0 = ffmin((_min[a] - r.origin()[a]) / r.direction()[a],
							(_max[a] - r.origin()[a]) / r.direction()[a]);
			float t1 = ffmax((_min[a] - r.origin()[a]) / r.direction()[a],
							(_max[a] - r.origin()[a]) / r.direction()[a]);
			tmin = ffmax(t0, tmin);
			tmax = ffmin(t1, tmax);
			if (tmax <= tmin)
				return false;
		}
		return true;
	}*/

	bool hit(const ray& r, float tmin, float tmax)const {
		for (int a = 0; a < 3; a++)
		{
			float invD = 1.0f / r.direction()[a];
			float t0 = (Min()[a] - r.origin()[a]) * invD;
			float t1 = (Max()[a] - r.origin()[a]) * invD;
			if (invD < 0.0f)
			{
				std::swap(t0, t1);
			}
			tmin = t0 > tmin ? t0 : tmin;
			tmax = t1 < tmax ? t1 : tmax;
			if (tmax <= tmin)
				return false;
		}
		return true;
	}
	vec3 _min;
	vec3 _max;
};

aabb surrounding_box(aabb box0, aabb box1) {
	vec3 _small(ffmin(box0.Min().x(), box1.Min().x()),
			   ffmin(box0.Min().y(), box1.Min().y()),
			   ffmin(box0.Min().z(), box1.Min().z()));
	vec3 _big  (ffmax(box0.Max().x(), box1.Max().x()),
			   ffmax(box0.Max().y(), box1.Max().y()),
			   ffmax(box0.Max().z(), box1.Max().z()));
	return aabb(_small, _big);
}
#endif // ! AABBH
