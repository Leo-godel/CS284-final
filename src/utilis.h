#pragma once
#include <ctime>
#include <cmath>
#include <cstdlib>
#include "Constants.h"
#include "Matrix2D.h"
#include "Vector2D.h"
#define MATRIX_EPSILON 1e-6

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
inline double bsplineSlope(double x) {
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

inline void svd(Matrix2D m, Matrix2D& U, Matrix2D& sig, Matrix2D& V){
    if (fabs(m(1, 0) - m(0, 1)) < MATRIX_EPSILON && fabs(m(1, 0)) < MATRIX_EPSILON) {
        U = Matrix2D(m(0, 0) < 0 ? -1 : 1, 0, 0, m(1, 1) < 0 ? -1 : 1);
        sig(0, 0) = fabs(m(0, 0)), sig(1, 1) = fabs(m(1, 1));
        V = Matrix2D();
    }
    else {
        double j = m(0, 0) * m(0, 0) + m(1, 0) * m(1, 0);
        double k = m(0, 1) * m(0, 1) + m(1, 1) * m(1, 1);
        double v_c = m(0, 0) * m(0, 1) + m(1, 0) * m(1, 1);

        if (fabs(v_c) < MATRIX_EPSILON) {
            double s1 = sqrt(j);
            double s2 = fabs(j - k) < MATRIX_EPSILON ? s1 : sqrt(k);
            sig(0, 0) = s1, sig(1, 1) = s2;
            V = Matrix2D();
            U = Matrix2D(m(0, 0) / s1, m(0, 1) / s2, m(1, 0) / s1, m(1, 1) / s2);
        }
        else {
            double jmk = j - k,
                    jpk = j + k,
                    root = sqrt(jmk * jmk + 4 * v_c * v_c),
                    eig = (jpk + root) / 2,
                    s1 = sqrt(eig),
                    s2 = fabs(root) < MATRIX_EPSILON ? s1 : sqrt((jpk - root) / 2);
            sig(0, 0) = s1, sig(1, 1) = s2;
            double v_s = eig - j,
                    len = sqrt(v_s * v_s + v_c * v_c);
            v_c /= len;
            v_s /= len;
            V = Matrix2D(v_c, -v_s, v_s, v_c);
            U = Matrix2D(
                    (m(0, 0) * v_c + m(0, 1) * v_s) / s1,
                    (m(0, 1) * v_c - m(0, 0) * v_s) / s2,
                    (m(1, 0) * v_c + m(1, 1) * v_s) / s1,
                    (m(1, 1) * v_c - m(1, 0) * v_s) / s2
            );
        }
    }
}