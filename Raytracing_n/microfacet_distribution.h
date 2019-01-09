#ifndef MICROFACET_DISTRIBUTION_H
#define MICROFACET_DISTRIBUTION_H
#include "vec3.h"
#include "mathf.h"
#include "reflection.h"
#include "geometry.h"

class MicrofacetDistribution {
public:
	MicrofacetDistribution();

	static vec3 BeckmannSample(const vec3& wi, float alpha_x, float alpha_y, float u1, float u2) {
		// 1. stretch wi
		vec3 wiStretched =
			unit_vector(vec3(alpha_x * wi.x(), alpha_y * wi.y(), wi.z()));

		// 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
		float slope_x, slope_y;
		BeckmannSample11(CosTheta(wiStretched), u1, u2, &slope_x, &slope_y);

		// rotate 
		float tmp = CosPhi(wiStretched) * slope_x - SinPhi(wiStretched) * slope_y;
		slope_y = SinPhi(wiStretched) * slope_x + CosPhi(wiStretched) * slope_y;
		slope_x = tmp;

		// 4. unstretch
		slope_x = alpha_x * slope_x;
		slope_y - alpha_y * slope_y;

		// 5. compute normal
		return unit_vector(vec3(-slope_x, -slope_y, 1.f));
	}

	static void BeckmannSample11(float cosThetaI, float u1, float u2, float* slope_x, float* slope_y) {
		/* Special case (normal incidence */
		if (cosThetaI > .9999) {
			float r = std::sqrt(-std::log(1.0f - u1));
			float sinPhi = sin(2 * M_PI * u2);
			float cosPhi = cos(2 * M_PI * u2);
			*slope_x = r * cosPhi;
			*slope_y = r * sinPhi;
			return;
		}

		/* The original inversion routine from the paper contained
		   discontunuties, which cases issues for QMC integration
		   and techniques like Kelemen-stype MLT. the following code
		   performs a numerical inversion with better behavior */
		float sinThetaI =
			sqrt(fmax((float)0, (float)1 - cosThetaI * cosThetaI));
		float tanThetaI = sinThetaI / cosThetaI;
		float cotThetaI = 1 / tanThetaI;

		/* Search interval -- everything is parameterized
		   is the Erf() domain */
		float a = -1, c = Erf(cosThetaI);
		float sample_x = fmax(u1, (float)1e-6f);

		/* Start with a good initial guess */
		// Float b = (1 - sample_x) * a * sample_x * c;

		/* We can do better (inverse of an Approximation computed in
		   Mathematica)*/
		float thetaI = acos(cosThetaI);
		float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
		float b = c * (1 + c) * pow(1 - sample_x, fit);

		/* Normalization factor for the CDF */
		static const float SQRT_PI_INV = 1.f / sqrt(M_PI);
		float normalization =
			1 /
			(1 + c + SQRT_PI_INV * tanThetaI * exp(-cotThetaI * cotThetaI));

		int it = 0;
		while (++it < 10)
		{
			/* Bisection criterion -- the oddly-looking
			   Boolean expression are intentional to check
			   for NaNs at little additional cost */
			if (!(b >= a && b <= c)) b = 0.5f * (a + c);

			/* Evaluate the CDF and its derivative
			   (i.e  the density function) */
			float invErf = ErfInv(b);
			float value =
				normalization *
				(1 + b + SQRT_PI_INV * tanThetaI* (exp(-invErf * invErf))) -
				sample_x;
			float derivative = normalization * (1 - invErf * tanThetaI);

			if (std::abs(value) < 1e-5f) break;

			/* Update bisection intervals */
			if (value > 0)
				c = b;
			else
				a = b;

			b -= value / derivative;
		}

		/* Now convert back into a slope value */
		*slope_x = ErfInv(b);

		/* Simulate Y component */
		*slope_y = ErfInv(2.0f * fmax(u2, (float)1e-6f) - 1.0f);
	}

	virtual float D(const vec3& wh) const = 0;
	virtual float Lambda(const vec3& w) const = 0;
	float G1(const vec3& w) const {
		return 1 / (1 + Lambda(w));
	}
	virtual float G(const vec3& wo, const vec3& wi) const {
		return 1 / (1 + Lambda(wo) + Lambda(wi));
	}

	virtual vec3 Sample_wh(const vec3& wo, const vec3& u) const = 0;

	float Pdf(const vec3& wo, const vec3& wh) const;

protected:
	MicrofacetDistribution(bool sampleVisibleArea) :
		sampleVisibleArea(sampleVisibleArea) {}

	const bool sampleVisibleArea;
};


float MicrofacetDistribution::Pdf(const vec3& wo, const vec3& wh) const {
	if (sampleVisibleArea)
		return D(wh) * G1(wo) * abs(dot(wo, wh)) / AbsCosTheta(wo);
	else
		return D(wh) * AbsCosTheta(wh);
}

class BeckmannDistribution : public MicrofacetDistribution {
public:
	static float RoughnessToAlpha(float roughness) {
		roughness = std::fmax(roughness, (float)1e-3);
		float x = std::log(roughness);
		return 1.62162f + 0.819955f *x + 0.1734f *x * x +
			0.0171201f *x * x * x + 0.000640711f *x *x *x*x;
	}

	BeckmannDistribution(float alphax, float alphay, bool samplevis = true)
		: MicrofacetDistribution(samplevis), alphax(alphax), alphay(alphay) {}
	float D(const vec3& wh) const;
	vec3 Sample_wh(const vec3& wo, const vec3& u) const;
private:
	float Lambda(const vec3& w) const;
	float alphax, alphay;
};

float BeckmannDistribution::D(const vec3& wh) const {
	float tan2Theta = Tan2Theta(wh);
	if (isinf(tan2Theta)) return 0.;
	float cos4Theta = Cos2Theta(wh) * Cos2Theta(wh);
	return exp(-tan2Theta * (Cos2Phi(wh) / (alphax * alphax) +
		Sin2Phi(wh) / (alphay * alphay))) /
		(M_PI * alphax * alphay * cos4Theta);
}

float BeckmannDistribution::Lambda(const vec3& w) const {
	float absTanTheta = abs(TanTheta(w));
	if (isinf(absTanTheta)) return 0;
	// Compute _alpha_ for direction _w_
	float alpha =
		sqrt(Cos2Phi(w) * alphax * alphax + Sin2Phi(w) * alphay * alphay);
	float a = 1 / (alpha * absTanTheta);
	if (a > 1.6f) return 0;
	return (1 - 1.259f * a + 0.396f *a * a) / (3.535f *a + 2.181f * a * a);
}

vec3 BeckmannDistribution::Sample_wh(const vec3 &wo, const vec3& u) const {
	if (!sampleVisibleArea)
	{
		float tan2theta, phi;
		if (alphax == alphay) {
			float logSample = std::log(1 - u[0]);
			tan2theta = -alphax * alphax * logSample;
			phi = u[1] * 2 * M_PI;
		}
		else {
			// Compute _tan2Theta_ and _phi_ for anisotropic Beckmann
			// distribution
			float logSample = std::log(1 - u[0]);
			phi = std::atan(alphay / alphax *
				std::tan(2 * M_PI * u[1] + 0.5f * M_PI));
			if (u[1] > 0.5f) phi += M_PI;
			float sinPhi = sin(phi), cosPhi = cos(phi);
			float alphax2 = alphax * alphax, alphay2 = alphay * alphay;
			tan2theta = -logSample / (cosPhi * cosPhi / alphax2 + sinPhi * sinPhi / alphay2);
		}

		// Mapsampled Beckmann angles to normal direction _wh_
		float cosTheta = 1 / sqrt(1 + tan2theta);
		float sinTheta = sqrt(fmaxf((float)0, 1 - cosTheta * cosTheta));
		vec3 wh = SphericalDirection(sinTheta, cosTheta, phi);
		if (!SameHemisphere(wo, wh)) wh = -wh;
		return wh;
	}
	else
	{
		vec3 wh;
		bool flip = wo.z() < 0;
		wh = BeckmannSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
		if (flip) wh = -wh;
		return wh;
	}
}


class GGXDistribution : public MicrofacetDistribution {
public :

	GGXDistribution() {}

};



#endif