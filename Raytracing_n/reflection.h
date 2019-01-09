#ifndef REFLECTION_H
#define REFLECTION_H

#include "vec3.h"
#include "common.h"

// BSDF Inline Function
inline float CosTheta(const vec3& w) { return w.z(); }
inline float Cos2Theta(const vec3 &w) { return w.z() * w.z(); }
inline float AbsCosTheta(const vec3 &w) { return std::abs(w.z()); }
inline float Sin2Theta(const vec3& w) {
	return std::fmaxf((float)0, (float)1 - Cos2Theta(w));
}

inline float SinTheta(const vec3& w) { return sqrt(Sin2Theta(w)); }
inline float TanTheta(const vec3& w) { return SinTheta(w) / CosTheta(w); }

inline float Tan2Theta(const vec3& w) { return Sin2Theta(w) / Cos2Theta(w); }

inline float CosPhi(const vec3& w) {
	float sinTheta = SinTheta(w);
	return (sinTheta == 0) ? 1 : Clamp(w.x() / sinTheta, -1, 1);
}

inline float SinPhi(const vec3& w) {
	float sinTheta = SinTheta(w);
	return (sinTheta == 0) ? 0 : Clamp(w.y() / sinTheta, -1, 1);
}

inline float Cos2Phi(const vec3& w) { return CosPhi(w) * CosPhi(w); }

inline float Sin2Phi(const vec3& w) { return SinPhi(w) * SinPhi(w); }

inline vec3 Reflect(const vec3& wo, const vec3& n) {
	return -wo + 2 * dot(wo, n) * n;
}

inline bool Refract(const vec3& wi, const vec3& n, float eta, vec3* wt) {
	// Compute $\cos \ theta_ \roman{t}$ using Snell's law
	float cosThetaI = dot(n, wi);
	float sin2ThetaI = fmax(float(0), float(1 - cosThetaI * cosThetaI));
	float sin2ThetaT = eta * eta * sin2ThetaI;

	// Handle total internal reflection for transmission
	if (sin2ThetaT >= 1) return false;
	float cosThetaT = sqrt(1 - sin2ThetaT);
	*wt = eta * -wi + (eta * cosThetaI - cosThetaT) * n;
	return true;
}

inline bool SameHemisphere(const vec3 &w, const vec3 &wp) {
	return w.z() * wp.z() > 0;
}

#endif
