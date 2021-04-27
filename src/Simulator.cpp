#include "Simulator.h"

void Simulator::update() {
    // P2G step
	grid->initGridMassVel();
//	grid->initGridVel();

	grid->computeForce();
	grid->updateVelocity();

	grid->updateDeformation();
	grid->updateParticlesVelocity();

	// update Particles
	scene->update();
}