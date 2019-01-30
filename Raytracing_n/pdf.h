#ifndef PDFH
#define PDFH
#include "vec3.h"
#include "onb.h"
#include "mathf.h"
#include "brdf.h"
#include "microfacet_distribution.h"
#include "rng.h"

inline vec3 random_cosine_direction() {
	float r1 = drand48();
	float r2 = drand48();
	float phi = 2 * M_PI*r1;
	float z = sqrt(1 - r2);
	float x = cos(phi) * 2 * sqrt(r2);
	float y = sin(phi) * 2 * sqrt(r2);
	return vec3(x, y, z);
}

static rng s_rng = rng();

class pdf {
public:
	virtual float value(const vec3& wo, const vec3& wi) const = 0;
	virtual vec3 generate(const vec3 &wo) const = 0;
	virtual ~pdf(){
	}
};

class cosine_pdf : public pdf {
public:
	cosine_pdf(const vec3& w): n(w) { uvw.build_from_w(w); }
	virtual float value(const vec3& wo, const vec3& wi) const {
		/*float cosine = dot(unit_vector(wi), uvw.w());
		if (cosine > 0)
			return cosine / M_PI;
		else
			return 0;*/

		float cosineO = dot(unit_vector(wo), n);
		float cosineI = dot(unit_vector(wi), n);
		if (cosineI * cosineO < 0)
			return std::abs(cosineI) / M_PI;
		else
			return 0;
	}
	virtual vec3 generate(const vec3 &wo) const {
		vec3 gdir = random_cosine_direction();
		if (dot(wo, n) > 0)
		{
			gdir.e[2] *= -1;
		}
		vec3 wi = uvw.local(gdir);

		return wi;
	}
	onb uvw;
	const vec3 n;
};

class onrennayar_pdf :public pdf {
public:
	onrennayar_pdf(const vec3& w) : n(w) { uvw.build_from_w(w); }
	virtual float value(const vec3& wo, const vec3& wi) const {
		//return 1;
		//vec3 half_vector = (unit_vector(wi) + unit_vector(wo)) / 2;
		float cosineO = dot(unit_vector(wo), n);
		float cosineI = dot(unit_vector(wi), n);
		if (cosineI * cosineO < 0)
			return std::abs(cosineI) / M_PI;
		else
			return 0;
	}

	virtual vec3 generate(const vec3& wo) const {
		vec3 gdir = random_cosine_direction();
		if (dot(wo, n) > 0)
		{
			gdir.e[2] *= -1;
		}
		vec3 wi = uvw.local(gdir);
		
		return wi;
	}
	onb uvw;
	const vec3 n;
};


class beckmann_pdf : public pdf {
public:
	beckmann_pdf(const BeckmannDistribution *distribution, const vec3& n) : distribution(distribution) {
		pdf_value = (float*)malloc(sizeof(float));
		uvw.build_from_w(n);
	}

	virtual ~beckmann_pdf()
	{
		//cout << "free pdf_value" << endl;
		free(pdf_value);
		//pdf_value = nullptr;
	}

	virtual float value(const vec3& wo, const vec3& wi) const { 

		return *pdf_value;
		//vec3 wwo = vec3(
		//	dot(unit_vector(wo), uvw.u()),
		//	dot(unit_vector(wo), uvw.v()),
		//	dot(unit_vector(wo), uvw.w()));
		//vec3 wwi = vec3(
		//	dot(unit_vector(wi), uvw.u()),
		//	dot(unit_vector(wi), uvw.v()),
		//	dot(unit_vector(wi), uvw.w())
		//);
		//vec3 wh = unit_vector(wwi + wwo);
		///*vec3 wo = unit_vector(r_in.direction());
		//vec3 wi = unit_vector(scattered.direction());*/

		//float cosThetaI = AbsCosTheta(wwi);
		//float cosThetaO = AbsCosTheta(wwo);

		//float cosine = dot(wh, wwi);
		////float F = schlick(cosine, 0.9);
		//float F = 0.8;
		//return distribution->D(wh) * distribution->G(wwo, wwi) * F / (4 * cosThetaI * cosThetaO);
	}
	virtual vec3 generate(const vec3 &wo) const {
		float u1 = s_rng.UniformFloat();
		float u2 = s_rng.UniformFloat();
		vec3 wwo = unit_vector(vec3(dot(-wo, uvw.u()), dot(-wo, uvw.v()), dot(-wo, uvw.w())));
		vec3 u(u1, u2, 0);
		vec3 wh = distribution->Sample_wh(wwo, u);
		vec3 wi = Reflect(unit_vector(wwo), wh);
		vec3 wwi = unit_vector(wi.x() * uvw.u() + wi.y() * uvw.v()+ wi.z() * uvw.w());
		/*if (SameHemisphere(wi, wwo))
		{
			*pdf_value = 0;
			return vec3(0.0);
		}*/

		
		
		*pdf_value = distribution->Pdf(wwo, wh) / (4 * dot(wwo, wh));
		if (SameHemisphere(wi, wwo))
		{
			*pdf_value = 0;
		}
		//wwi *= distribution->Pdf(unit_vector(wwo), wh) / (4 * dot(unit_vector(wwo), wh));
		return wwi;
	}
	float *pdf_value;
	const BeckmannDistribution *distribution;
	onb uvw;
};


class hitable_pdf :public pdf {
public:
	hitable_pdf(hitable *p, const vec3& origin) : ptr(p), o(origin) {}
	virtual float value(const vec3& wo, const vec3& wi) const {
		return ptr->pdf_value(o, wi);
	}
	virtual vec3 generate(const vec3 &wo) const {
		return ptr->random(o);
	}
	~hitable_pdf() {}
	vec3 o;
	hitable *ptr;
};

class mixture_pdf : public pdf {
public:
	mixture_pdf(pdf *p0, pdf* p1) { p[0] = p0; p[1] = p1; rand = drand48(); }
	virtual float value(const vec3& wo, const vec3& wi) const {
	/*	if (rand < 0.5)
			return p[0]->value(direction);
		else
			return p[1]->value(direction);*/
		return 0.5 * p[0]->value(wo, wi) + 0.5 * p[1]->value(wo, wi);
	}
	virtual vec3 generate(const vec3 &wo) const {
		if (drand48() < 0.5)
		//if (rand < 0.5)
			return p[0]->generate(wo);
		else
			return p[1]->generate(wo);
	}

	pdf* p[2];
	float rand;
};
#endif