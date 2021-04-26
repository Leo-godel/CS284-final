#pragma once
#include <ctime>
#include <cstdlib>

inline double randomLR(double left, double right) {
	return left + rand() / (double)(RAND_MAX / (right - left));
}
