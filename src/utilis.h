#pragma once
#include <ctime>
#include <cstdlib>
#include "Constants.h"
#include "Matrix2D.h"
#include "Vector2D.h"

inline double randomLR(double left, double right) {
	return left + rand() / (double)(RAND_MAX / (right - left));
}

inline double bspline(double x) {
	x = fabs(x);
	double w;
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