#include "Simulator.h"

void Simulator::update() {
    // P2G step
	// cout << "0" << endl;
	grid->G2P2G();
	// grid->initGridMassVel();
	// cout << "1" << endl;
//	grid->initGridVel();

	grid->computeForce();
	grid->updateVelocity();

	// grid->updateDeformation();
	// grid->updateParticlesVelocity();
	// cout << "2" << endl;
	// grid->updateParticlesPosition();
	// cout << "3" << endl;
}