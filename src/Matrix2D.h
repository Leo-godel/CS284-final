#pragma once

#include "Vector2D.h"
#include <cmath>

class Matrix2D
{
public:

	// Default Matrix = Identity Matrix
	Matrix2D(double val = 1.) {
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				(*this)(i, j) = (i == j) ? val : 0.;
	}
	Matrix2D(double* data) {
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				(*this)(i, j) = data[i * 2 + j];
	}
	Matrix2D(double m00, double m01,
		double m10, double m11) {
		(*this)(0, 0) = m00; (*this)(0, 1) = m01;
		(*this)(1, 0) = m10; (*this)(1, 1) = m11;
	}

	// set all values to val (detaul: 0.)
	void zero(double val = 0.);

	// determinant
	double det(void) const;

	// Return the column vector
	Vector2D& column(int i);
	const Vector2D& column(int i) const;

	// Transpose
	Matrix2D T(void) const;

	// Inverse
	Matrix2D inv(void) const;

	// accesses element (i,j) of A using 0-based indexing
	double& operator()(int i, int j);
	const double& operator()(int i, int j) const;

	// accesses the ith column of A
	Vector2D& operator[](int i);
	const Vector2D& operator[](int i) const;

	// increments by B
	void operator+=(const Matrix2D& B);

	// A+c
	Matrix2D operator+(const Matrix2D& B) const;

	// returns -A
	Matrix2D operator-(void) const;

	// returns A-B
	Matrix2D operator-(const Matrix2D& B) const;

	// returns c*A
	Matrix2D operator*(double c) const;

	// returns A*B
	Matrix2D operator*(const Matrix2D& B) const;

	// returns A*x
	Vector2D operator*(const Vector2D& x) const;

	// divides each element by x
	void operator/=(double x);

protected:
	Vector2D entries[2];
};

