#ifndef PARTICLE_H
#define PARTICLE_H


#include <cstring>
#include "Vector2D.h"
#include "Matrix2D.h"
#include <cuda.h>
#include <cuda_runtime.h>


class Particle
{
public:
	double mass, density, volume;
	Vector2D pos, vel;
	Matrix2D velocity_gradient;
	Matrix2D deformation_gradient, plastic_deformation, elastic_deformation;
	Matrix2D stress;
	
	// assigned grid index
	Vector2D grid_p;
	// weight value for nearest 16 grid nodes
	double weights[16];
	Vector2D weight_gradient[16];

	// TODO: need more variables here like deformation;

	__host__ __device__ Particle() {};
	
	__host__ __device__ Particle(const Vector2D& pos, const Vector2D& vel, double mass) :
		pos(pos), vel(vel), mass(mass) {
		// TODO: more varibles means re-write init function
		this->deformation_gradient = Matrix2D();  // init as identity matrix
        this->plastic_deformation = Matrix2D();  // init as identity matrix
        this->elastic_deformation = Matrix2D();  // init as identity matrix
		this->volume = 1;  // init volume is 1 (a unit volume)
		memset(weights, 0, sizeof(double) * 16);
        memset(weight_gradient, 0, sizeof(Vector2D) * 16);
	}
};

#endif