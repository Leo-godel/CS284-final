#ifndef MATRIX_2D_H
#define MATRIX_2D_H

#include "Vector2D.h"
#include <cmath>
#include <cuda.h>
#include <cuda_runtime.h>

class Matrix2D
{
public:

	// Default Matrix = Identity Matrix
	__host__ __device__ Matrix2D(double val = 1.) {
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				(*this)(i, j) = (i == j) ? val : 0.;
	}
	__host__ __device__ Matrix2D(double* data) {
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				(*this)(i, j) = data[i * 2 + j];
	}
	__host__ __device__ Matrix2D(double m00, double m01,
		double m10, double m11) {
		(*this)(0, 0) = m00; (*this)(0, 1) = m01;
		(*this)(1, 0) = m10; (*this)(1, 1) = m11;
	}

	// set all values to val (detaul: 0.)
	__host__ __device__ void zero(double val = 0.);

	// determinant
	__host__ __device__ double det(void) const;

	// Return the column vector
	__host__ __device__ Vector2D& column(int i);
	__host__ __device__ const Vector2D& column(int i) const;

	// Transpose
	__host__ __device__ Matrix2D T(void) const;

	// Inverse
	__host__ __device__ Matrix2D inv(void) const;

	// accesses element (i,j) of A using 0-based indexing
	__host__ __device__ double& operator()(int i, int j);
	__host__ __device__ const double& operator()(int i, int j) const;

	// accesses the ith column of A
	__host__ __device__ Vector2D& operator[](int i);
	__host__ __device__ const Vector2D& operator[](int i) const;

	// increments by B
	__host__ __device__ void operator+=(const Matrix2D& B);

	// A+c
	__host__ __device__ Matrix2D operator+(const Matrix2D& B) const;

	// returns -A
	__host__ __device__ Matrix2D operator-(void) const;

	// returns A-B
	__host__ __device__ Matrix2D operator-(const Matrix2D& B) const;

	// returns c*A
	__host__ __device__ Matrix2D operator*(double c) const;

	// returns A*B
	__host__ __device__ Matrix2D operator*(const Matrix2D& B) const;

	// returns A*x
	__host__ __device__ Vector2D operator*(const Vector2D& x) const;

	// divides each element by x
	__host__ __device__ void operator/=(double x);

protected:
	Vector2D entries[2];
};

#endif