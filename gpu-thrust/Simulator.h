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

		//cout << grid->node_size.x << " " << grid->node_size.y << endl;
		//cout << grid->size.x << " " << grid->size.y << endl;

		grid->initGridMassVel();
		grid->initVolumes();
	}

	void update();
};

#endif