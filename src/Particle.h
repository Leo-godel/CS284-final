#pragma once

#include "Vector2D.h"
#include "Matrix2D.h"

class Particle
{
public:
	double mass, density, volume;
	Vector2D pos, vel;
	// Lame parameters
	double lambda, mu;
	// assigned grip position
	Vector2D grid_p;
	// weight value for nearest 16 grid nodes
	double weights[16];
	Vector2D weight_gradient[16];

	// TODO: need more variables here like deformation;

	Particle() {};
	
	Particle(const Vector2D& pos, const Vector2D& vel, double mass, double lambda, double mu) :
		pos(pos), vel(vel), mass(mass), lambda(lambda), mu(mu) {
		// TODO: more varibles means re-write init function
	}

	// TODO: position update
	void updatePosition();
	// TODO: inner force update
	void updateForce();
	// TODO: calculate stree tensor
	void calcStress();
};

