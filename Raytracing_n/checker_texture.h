#pragma once
#ifndef CHECKER_TEXTUREH
#define CHCKER_TEXTUREH
#include "texture.h"

class checker_texture : public texture {
public:
	checker_texture() {}
	checker_texture(texture *t0, texture*t1) :even(t0), odd(t1) {}
	virtual vec3 value(float u, float v, const vec3& p)const {
		float sines = sin(10 * p.x()) * sin(10 * p.y())* sin(10 * p.z());
		if (sines < 0)
			return odd->value(u, v, p);
		else
			return even->value(u, v, p);
	}

	texture *odd;
	texture *even;
};

#endif