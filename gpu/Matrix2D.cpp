#include "Matrix2D.h"

double& Matrix2D::operator() (int i, int j) {
	return entries[j][i];
}
const double& Matrix2D::operator() (int i, int j) const {
	return entries[j][i];
}

Vector2D& Matrix2D::operator[] (int j) {
	return entries[j];
}
const Vector2D& Matrix2D::operator[] (int j) const {
	return entries[j];
}

void Matrix2D::zero(double val) {
	entries[0] = entries[1] = Vector2D(val, val);
}

double Matrix2D::det(void) const {
	const Matrix2D& A(*this);

	return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
}

Matrix2D Matrix2D::T(void) const {
	const Matrix2D& A(*this);
	Matrix2D B;
	for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j)
            B(i, j) = A(j, i);
	}
	return B;
}

Matrix2D Matrix2D::inv(void) const {
	const Matrix2D& A(*this);
	double d = det();
	return Matrix2D(A(1, 1) / d, -A(1, 0) / d, -A(0, 1) / d, A(0, 0) / d);
}


Vector2D& Matrix2D::column(int i) {
	return entries[i];
}

Matrix2D Matrix2D::operator-(void) const {
	const Matrix2D& A(*this);
	Matrix2D B;
	B[0] = -A[0];
	B[1] = -A[1];
	
	return B;
}

Matrix2D Matrix2D::operator-(const Matrix2D& B) const {
	const Matrix2D& A(*this);
	Matrix2D C;
	C[0] = A[0] - B[0];
	C[1] = A[1] - B[1];

	return C;
}

void Matrix2D::operator+=(const Matrix2D& B) {
	Matrix2D& A(*this);

	A[0] += B[0];
	A[1] += B[1];
}

Matrix2D Matrix2D::operator*(double c) const {
	const Matrix2D& A(*this);
	Matrix2D B;

	B[0] = A[0] * c;
	B[1] = A[1] * c;
	return B;
}

Matrix2D Matrix2D::operator*(const Matrix2D& B) const {
	const Matrix2D& A(*this);
	Matrix2D C;

	C[0] = A * B[0];
	C[1] = A * B[1];

	return C;
}

Vector2D Matrix2D::operator*(const Vector2D& x) const {
	return x.x * entries[0] + x.y * entries[1];
}

void Matrix2D::operator/=(double x) {
	Matrix2D& A(*this);
	double rx = 1. / x;

	A[0] *= rx;
	A[1] *= rx;
}

Matrix2D Matrix2D::operator+(const Matrix2D& B) const {
	const Matrix2D& A(*this);
	Matrix2D C;
	C[0] = A[0] + B[0];
	C[1] = A[1] + B[1];

	return C;
}