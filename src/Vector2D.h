#pragma once

#include <cmath>

class Vector2D
{
public:
	double x, y;

	// Constructor
	Vector2D() : x(0.), y(0.) {}
	Vector2D(double x, double y) : x(x), y(y) {}
	Vector2D(const Vector2D& v) : x(v.x), y(v.y) {}

	//Operators
	// returns reference to the specified component (0-based indexing: x, y)
	inline double& operator[] (const int& index) {
		return (&x)[index];
	}

	// returns const reference to the specified component (0-based indexing: x, y)
	inline const double& operator[] (const int& index) const {
		return (&x)[index];
	}

	// additive inverse
	inline Vector2D operator-(void) const {
		return Vector2D(-x, -y);
	}

	// addition
	inline Vector2D operator+(const Vector2D& v) const {
		Vector2D u = *this;
		u += v;
		return u;
	}

	// subtraction
	inline Vector2D operator-(const Vector2D& v) const {
		Vector2D u = *this;
		u -= v;
		return u;
	}

	// right scalar multiplication
	inline Vector2D operator*(double r) const {
		Vector2D vr = *this;
		vr *= r;
		return vr;
	}
	inline Vector2D operator*(const Vector2D& v) const {
		Vector2D vr = *this;
		vr.x *= v.x;
		vr.y *= v.y;
		return vr;
	}


	// scalar division
	inline Vector2D operator/(double r) const {
		Vector2D vr = *this;
		vr /= r;
		return vr;
	}

	inline Vector2D operator/(Vector2D v) const {
		Vector2D vr = *this;
		vr.x = vr.x / v.x;
		vr.y = vr.y / v.y;
		return vr;
	}

	// add v
	inline void operator+=(const Vector2D& v) {
		x += v.x;
		y += v.y;
	}

	// subtract v
	inline void operator-=(const Vector2D& v) {
		x -= v.x;
		y -= v.y;
	}

	// scalar multiply by r
	inline void operator*=(double r) {
		x *= r;
		y *= r;
	}

	// scalar divide by r
	inline void operator/=(double r) {
		x /= r;
		y /= r;
	}
	inline void operator/=(Vector2D v) {
		*this = *this / v;
	}

	/**
	 * Returns norm.
	 */
	inline double norm(void) const {
		return sqrt(x * x + y * y);
	}

	/**
	 * Returns norm squared.
	 */
	inline double norm2(void) const {
		return x * x + y * y;
	}

	/**
	 * Returns unit vector parallel to this one.
	 */
	inline Vector2D unit(void) const {
		return *this / this->norm();
	}


}; // clasd Vector2D

// left scalar multiplication
inline Vector2D operator*(double r, const Vector2D& v) {
	return v * r;
}

// inner product
inline double dot(const Vector2D& v1, const Vector2D& v2) {
	return v1.x * v2.x + v1.y * v2.y;
}

// cross product
inline double cross(const Vector2D& v1, const Vector2D& v2) {
	return v1.x * v2.y - v1.y * v2.x;
}


