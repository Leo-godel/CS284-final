#include "Simulator.h"

void Simulator::update() {
    // P2G step
	grid->initGridMassVel();
//	grid->initGridVel();

	grid->computeForceVel();

	grid->updateParticles();

	// update Particles
	scene->update();

    // G2P step
//    grid->updateVelocities();
//
//    //Update particle data
//    snow->update();
}