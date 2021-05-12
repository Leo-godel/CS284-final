#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>

class Vector2D
{
public:
	double x, y;

	// Constructor
	__host__ __device__ Vector2D() : x(0.), y(0.) {}
	__host__ __device__ Vector2D(double x, double y) : x(x), y(y) {}
	__host__ __device__ Vector2D(const Vector2D& v) : x(v.x), y(v.y) {}

	//Operators
	// returns reference to the specified component (0-based indexing: x, y)
	__host__ __device__ inline double& operator[] (const int& index) {
		return (&x)[index];
	}

	// returns const reference to the specified component (0-based indexing: x, y)
	__host__ __device__ inline const double& operator[] (const int& index) const {
		return (&x)[index];
	}

	// additive inverse
	__host__ __device__ inline Vector2D operator-(void) const {
		return Vector2D(-x, -y);
	}

	// addition
	__host__ __device__  inline Vector2D operator+(const Vector2D& v) const {
		Vector2D u = *this;
		u += v;
		return u;
	}

	// subtraction
	__host__ __device__  inline Vector2D operator-(const Vector2D& v) const {
		Vector2D u = *this;
		u -= v;
		return u;
	}

	// right scalar multiplication
	__host__ __device__  inline Vector2D operator*(double r) const {
		Vector2D vr = *this;
		vr *= r;
		return vr;
	}
	__host__ __device__  inline Vector2D operator*(const Vector2D& v) const {
		Vector2D vr = *this;
		vr.x *= v.x;
		vr.y *= v.y;
		return vr;
	}


	// scalar division
	__host__ __device__  inline Vector2D operator/(double r) const {
		Vector2D vr = *this;
		vr /= r;
		return vr;
	}

	__host__ __device__  inline Vector2D operator/(Vector2D v) const {
		Vector2D vr = *this;
		vr.x = vr.x / v.x;
		vr.y = vr.y / v.y;
		return vr;
	}

	// add v
	__host__ __device__  inline void operator+=(const Vector2D& v) {
		x += v.x;
		y += v.y;
	}

	// subtract v
	__host__ __device__  inline void operator-=(const Vector2D& v) {
		x -= v.x;
		y -= v.y;
	}

	// scalar multiply by r
	__host__ __device__  inline void operator*=(double r) {
		x *= r;
		y *= r;
	}

	// scalar divide by r
	__host__ __device__  inline void operator/=(double r) {
		x /= r;
		y /= r;
	}
	__host__ __device__  inline void operator/=(Vector2D v) {
		*this = *this / v;
	}

	/**
	 * Returns norm.
	 */
	__host__ __device__  inline double norm(void) const {
		return sqrt(x * x + y * y);
	}

	/**
	 * Returns norm squared.
	 */
	__host__ __device__  inline double norm2(void) const {
		return x * x + y * y;
	}

	/**
	 * Returns unit vector parallel to this one.
	 */
	__host__ __device__  inline Vector2D unit(void) const {
		return *this / this->norm();
	}


}; // clasd Vector2D

// left scalar multiplication
__host__ __device__ inline Vector2D operator*(double r, const Vector2D& v) {
	return v * r;
}

// inner product
__host__ __device__ inline double dot(const Vector2D& v1, const Vector2D& v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

// cross product
__host__ __device__ inline double cross(const Vector2D& v1, const Vector2D& v2) {
	return v1.x * v2.y - v1.y * v2.x;
}


#endif