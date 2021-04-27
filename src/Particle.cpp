#include "Particle.h"
#include "Constants.h"

void Particle::updateForce() {

}

void Particle::updatePosition() {
	this->pos += TIMESTEP * this->vel;
}

Matrix2D Particle::calcStress() {
	return Matrix2D();
}