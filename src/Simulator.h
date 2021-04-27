#pragma once

#include "Grid.h"
#include "Scene.h"

class Simulator
{
public:
	Grid* grid;
	Scene* scene;

	Simulator() {}
	Simulator(Scene* s) :
		scene(s) {
		grid = new Grid(Vector2D(0, 0), Vector2D(WIN_METERS_X, WIN_METERS_Y), Vector2D(256, 128), s);

		grid->initGridMassVel();
		grid->initVolumes();
	}

	void update();
};

