#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "Grid.h"
#include "Scene.h"
#include <cuda.h>
#include <cuda_runtime.h>


class Simulator
{
public:
	Grid* grid;
	Scene* scene;

	__host__ Simulator() {}
	__host__ Simulator(Scene* s) :
		scene(s) {
		grid = new Grid(Vector2D(0, 0), Vector2D(WIN_METERS_X, WIN_METERS_Y), Vector2D(256, 128), s);

		// cout << "0" << endl;
		grid->initGridMassVel();
		grid->initVolumes();
		// cout << "1" << endl;
		grid->computeForce();
		grid->updateVelocity();
	}

	void update();
};

#endif