#ifndef MATERIAL_H
#define MATERIAL_H

#include "common.h"
#include "reflection.h"
#include "ray.h"
#include "hitable.h"
#include "texture.h"
#include "onb.h"
#include "brdf.h"
#include "pdf.h"
#include "microfacet_distribution.h"

float schlick(float cosine, float ref_idx)
{
	float r0 = (1 - ref_idx) / (1 + ref_idx);
	r0 = r0 * r0;
	return r0 + (1 - r0)*pow((1 - cosine), 5);
}

bool refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted) {
	vec3 uv = unit_vector(v);
	float dt = dot(uv, n);
	float discriminant = 1.0 - ni_over_nt * ni_over_nt *(1 - dt * dt);
	if (discriminant > 0)
	{
		refracted = ni_over_nt * (uv - n * dt) - n * sqrt(discriminant);
		return true;
	}
	else
		return false;
}

vec3 reflect(const vec3& v, const vec3& n) {
	return v - 2 * dot(v, n)*n;
}
float rand01()
{
	//return (rand() % 100) / (float)100;
	return drand48();
}

vec3 random_in_unit_sphere()
{
	vec3 p;
	do {
		p = 2.0 * vec3(rand01(), rand01(), rand01()) - vec3(1.0, 1.0, 1.0);
	} while (dot(p,p) >= 1.0);
	return p;
}

vec3 random_on_unit_sphere()
{
	vec3 p;
	do {
		p = 2.0*vec3(drand48(), drand48(), drand48()) - vec3(1.0, 1.0, 1.0);
	} while (dot(p, p) >= 1.0);
	return unit_vector(p);
}




struct scatter_record
{
	ray specular_ray;
	bool is_specular;
	vec3 attenuation;
	pdf *pdf_ptr;
};

enum BXDFType
{
	BSDF_REFLECTION = 1 << 0,
	BSDF_TRANSMISSION = 1 << 1,
	BSDF_DIFFUSE = 1 << 2,
	BSDF_GLOSSY = 1 << 3,
	BSDF_SPECULAR = 1 << 4,
	BSDF_ALL = BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR,
};

class material
{
public:
	virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec) const = 0;
	virtual float scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {
		return false;
	}
	virtual vec3 emitted(const ray& r_in, const hit_record& rec, float u, float v, const vec3& p)
	{
		return vec3(0, 0, 0);
	}
};

class lambertian : public material {
public:
	lambertian(texture *a) : albedo(a) {}

	
	float scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {
		float consine = dot(rec.normal, unit_vector(scattered.direction()));
		if (consine < 0) consine = 0;
		return consine / M_PI;
		//return 1 / M_PI;
	}
	virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec) const {
		srec.is_specular = false;
		srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
		srec.pdf_ptr = new cosine_pdf(hrec.normal);
		return true;
	}

	texture *albedo;
};

class disney_diffuse :public material {
public:
	disney_diffuse(texture *a) :albedo(a) {}

	float scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {

	}

	texture *albedo;
};

class orennayar : public material {
public:
	orennayar(texture *a, float sigma) : albedo(a) {
		sigma = sigma / 180 * M_PI;
		A = 1 - 0.5 * sigma * sigma / (sigma * sigma + 0.33);
		B = 0.45 * sigma * sigma / (sigma * sigma + 0.09);
	}
	float scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {
		float consine = dot(rec.normal, unit_vector(scattered.direction()));
		if (consine < 0) consine = 0;
		return consine / M_PI;
	}

	virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec) const {
		srec.is_specular = false;
		srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
		srec.pdf_ptr = new onrennayar_pdf(hrec.normal, A, B);
		return true;
	}
	texture* albedo;
	float A, B;
};

class beckmann : public material {
public:
	beckmann(texture *a, float roughx, float roughy) : albedo(a), roughx(roughx), roughy(roughy) {
		float alphax = BeckmannDistribution::RoughnessToAlpha(roughx);
		float alphay = BeckmannDistribution::RoughnessToAlpha(roughy);
		distribution = new BeckmannDistribution(alphax, alphay, true);
	}


	float scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {
		onb uvw;
		uvw.build_from_w(rec.normal);

		
		vec3 wo = unit_vector(vec3(
			dot(unit_vector(-r_in.direction()), uvw.u()), 
			dot(unit_vector(-r_in.direction()), uvw.v()), 
			dot(unit_vector(-r_in.direction()), uvw.w())
		));
		vec3 wi = unit_vector(vec3(
			dot(unit_vector(scattered.direction()), uvw.u()), 
			dot(unit_vector(scattered.direction()), uvw.v()), 
			dot(unit_vector(scattered.direction()), uvw.w())
		));
		vec3 wh = unit_vector(wi + wo);
		/*vec3 wo = unit_vector(r_in.direction());
		vec3 wi = unit_vector(scattered.direction());*/

		float cosThetaI = AbsCosTheta(wi);
		float cosThetaO = AbsCosTheta(wo);

		float F = 1;
		//return distribution->D(wh) * distribution->G(wo, wi) * F / (4 * cosThetaI * cosThetaO);
		return distribution->Pdf(wo, wh) / (4 * dot(wo, wh));
	}

	virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec) const {
		srec.is_specular = false;
		srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
		srec.pdf_ptr = new beckmann_pdf(distribution, hrec.normal);
		return true;
	}

	

	texture* albedo;
	const float roughx, roughy;	// alphax 垂直于x轴方向 alphay 垂直于y轴方向
	BeckmannDistribution *distribution;
};

class brdfmaterial : public material {
public :
	brdfmaterial(const char* brdf_filename, const vec3& a): albedo(a) {
		brdf_reader = new brdf();
		brdf_reader->read_brdf(brdf_filename, brdf_data);
	}
	float scattering_pdf(const ray& r_in, const hit_record& hrec, const ray& scattered)const {
		
		float cosine = dot(hrec.normal, unit_vector(scattered.direction()));
		if (cosine < 0) cosine = 0;
		return cosine / M_PI;
		//return 1;
	}
	virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec) const {
		vec3 reflected = reflect(unit_vector(r_in.direction()), hrec.normal);
		srec.specular_ray = ray(hrec.p, reflected);
		float cosine = dot(-unit_vector(r_in.direction()), hrec.normal);
		float theta_in = acos(cosine);
		vec3 in_right = cross(-r_in.direction(), vec3(0, 1, 0));

		vec3 phi_vector = hrec.normal + unit_vector(r_in.direction()) * cosine;
		float phi_in = acos(dot(unit_vector(in_right), unit_vector(phi_vector)));
		float phi_out;
		if (phi_in < 0)
			phi_in += 2 * M_PI;
		if (phi_in < M_PI)
			phi_out = phi_in + M_PI;
		else
			phi_out = phi_in - M_PI;
		double red, green, blue;
		brdf_reader->lookup_brdf_val(brdf_data, theta_in, phi_in, theta_in, phi_in, red, green, blue);
		srec.attenuation = albedo;// vec3(red, green, blue) * M_PI * M_PI / cosine;
		srec.is_specular = false;
		//srec.pdf_ptr = new reflect_pdf(hrec.normal, r_in.direction());
		return true;
	}

	brdf* brdf_reader;
	double* brdf_data;
	vec3 albedo;
};

class metal : public material {
public:
	metal(const vec3& a, float f) : albedo(a) { if (f < 1) fuzz = f; else fuzz = 1; }
	

	virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record& srec)const {
		vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
		srec.specular_ray = ray(rec.p, reflected + fuzz * random_in_unit_sphere());

		srec.attenuation = albedo;
		srec.is_specular = true;
		srec.pdf_ptr = 0;
		return true;
	}
	vec3 albedo;
	brdf* brdfreader;
	double* brdf_data;
	float fuzz;
};
//
//class toon : public material {
//public:
//	toon(const texture* a) : albedo(a){}
//	   float scattering_pdf(const ray& r_in, const hit_record& rec, const ray& scattered) const {
//		   float consine = dot(rec.normal, unit_vector(scattered.direction()));
//		   if (consine < 0) consine = 0;
//		   return consine / M_PI;
//		   //return 1 / M_PI;
//	   }
//	   virtual bool scatter(const ray& r_in, const hit_record& hrec, scatter_record& srec) const {
//		   srec.is_specular = false;
//		   srec.attenuation = albedo->value(hrec.u, hrec.v, hrec.p);
//		   srec.pdf_ptr = new cosine_pdf(hrec.normal);
//		   return true;
//	   }
//
//	   texture *albedo;
//};

class dielectric : public material {
public:
	dielectric(float ri) : ref_idx(ri) {}		//refractive indices (typically air = 1, glass = 1.3-1.7, diamond = 2.4)
	virtual bool scatter(const ray&r_in, const hit_record& rec, scatter_record& srec)const {
		srec.is_specular = true;
		srec.pdf_ptr = 0;
		srec.attenuation = vec3(1.0, 1.0, 1.0);
		vec3 outward_normal;
		vec3 reflected = reflect(r_in.direction(), rec.normal);
		float ni_over_nt;
		vec3 refracted;
		float reflect_prob;
		float cosine;
		if (dot(r_in.direction(), rec.normal) > 0) {		// out
			outward_normal = -rec.normal;
			ni_over_nt = ref_idx;			// sint / sini

			cosine = dot(r_in.direction(), rec.normal) / r_in.direction().length();
		}
		else {												// in
			outward_normal = rec.normal;
			ni_over_nt = 1.0 / ref_idx;
			cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
		}
		if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted))
		{
			reflect_prob = schlick(cosine, ref_idx);
		}
		else {
			srec.specular_ray = ray(rec.p, reflected);
			reflect_prob = 1.0;
		}
		if (rand01() < reflect_prob)
		{
			srec.specular_ray = ray(rec.p, reflected);
		}
		else
		{
			srec.specular_ray = ray(rec.p, refracted);
			//srec.attenuation = vec3(1 - reflect_prob, 1 - reflect_prob, 1 - reflect_prob);
		}
		return true;
	}

	// fresnl reflection equation  
	// energy transmitted by a dielectric is 1 - Fr
	float frdielectric(float cosine, float etai, float etat) {
		if (cosine < -1)
			cosine = -1;
		if (cosine > 1)
			cosine = 1;
		float cosine_t = sqrt(ffmax(0.0, 1.0 - etai * etai * (1 - cosine * cosine) / (etat * etat)));
		float rparl = ((etat * cosine) - (etai * cosine_t)) / ((etat * cosine) + (etai * cosine_t));
		float rperp = ((etai * cosine) - (etat * cosine_t)) / ((etai * cosine) + (etat * cosine_t));
		return (rparl * rparl + rperp * rperp) / 2;
	}
	float ref_idx;
};

class diffuse_light :public material {
public:
	diffuse_light(texture *a) :emit(a) {}
	virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record &srec) const { return false; }
	virtual vec3 emitted(float u, float v, const vec3& p) const {
		return emit->value(u, v, p);
	}
	virtual vec3 emitted(const ray& r_in, const hit_record& rec, float u, float v, const vec3& p)
	{
		if (dot(rec.normal, r_in.direction()) < 0.0)
			return emit->value(u, v, p);
		else
			return vec3(0, 0, 0);
	}
	texture* emit;
};

// isotropic : 各向同性
class isotropic : public material {
public:
	isotropic(texture* a) : albedo(a) {}
	virtual bool scatter(const ray& r_in, const hit_record& rec, scatter_record &srec) const {
		srec.is_specular = true;
		srec.specular_ray = ray(rec.p, random_in_unit_sphere());
		srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
		return true;
	}
	texture *albedo;
};
#endif // MATERIAL_H
