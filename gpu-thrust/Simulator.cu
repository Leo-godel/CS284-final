#include "Simulator.h"

void Simulator::update() {
    // P2G step

	grid->initGridMassVel();

	grid->computeForce();
	grid->updateVelocity();

	grid->updateDeformation();
	grid->updateParticlesVelocity();
	grid->updateParticlesPosition();
}