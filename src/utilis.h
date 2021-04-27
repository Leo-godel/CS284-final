#pragma once
#include <ctime>
#include <cmath>
#include <cstdlib>
#include "Constants.h"
#include "Matrix2D.h"
#include "Vector2D.h"

inline double randomLR(double left, double right) {
	return left + rand() / (double)(RAND_MAX / (right - left));
}

inline double bspline(double x) {
	x = fabs(x);
	float w;
	if (x < 1)
		w = x * x * (x / 2 - 1) + 2 / 3.0;
	else if (x < 2)
		w = x * (x * (-x / 6 + 1) - 2) + 4 / 3.0;
	else return 0;
	//Clamp between 0 and 1... if needed
	if (w < BSPLINE_EPSILON) return 0;
	return w;
}
//Slope of interpolation function
inline double bsplineSlope(float x) {
	double abs_x = fabs(x), w;
	if (abs_x < 1)
		return 1.5 * x * abs_x - 2 * x;
	else if (x < 2)
		return -x * abs_x / 2 + 2 * x - 2 * x / abs_x;
	else return 0;
	//Clamp between -2/3 and 2/3... if needed
}

inline Matrix2D outer_product(Vector2D& a, Vector2D& b) {
	return Matrix2D(a.x * b.x, a.x * b.y, a.y * b.x, a.y * b.y);
}

inline void polarDecomp(Matrix2D m, Matrix2D& R, Matrix2D& S) {
	auto x = m(0, 0) + m(1, 1);
	auto y = m(1, 0) - m(0, 1);
	auto scale = 1.0 / sqrt(x * x + y * y);
	auto c = x * scale, s = y * scale;
	R(0, 0) = c;
	R(0, 1) = -s;
	R(1, 0) = s;
	R(1, 1) = c;
	S = R.T() * m;
}

inline void svd(Matrix2D m, Matrix2D& U, Matrix2D& sig, Matrix2D& V) {
	Matrix2D S;
	polarDecomp(m, U, S);
	double c, s;
	if (fabs(S(0, 1)) < 1e-6) {
		sig = S;
		c = 1;
		s = 0;
	}
	else {
		auto tao = 0.5 * (S(0, 0) - S(1, 1));
		auto w = sqrt(tao * tao + S(0, 1) * S(0, 1));
		auto t = tao > 0 ? S(0, 1) / (tao + w) : S(0, 1) / (tao - w);
		c = 1.0 / sqrt(t * t + 1);
		s = -t * c;
		sig(0, 0) = pow(c, 2) * S(0, 0) - 2 * c * s * S(0, 1) + pow(s, 2) * S(1, 1);
		sig(1, 1) = pow(s, 2) * S(0, 0) + 2 * c * s * S(0, 1) + pow(c, 2) * S(1, 1);
	}
	if (sig(0, 0) < sig(1, 1)) {
		std::swap(sig(0, 0), sig(1, 1));
		V(0, 0) = -s;
		V(0, 1) = -c;
		V(1, 0) = c;
		V(1, 1) = -s;
	}
	else {
		V(0, 0) = c;
		V(0, 1) = -s;
		V(1, 0) = s;
		V(1, 1) = c;
	}
	V = V.T();
	U = U * V;
}