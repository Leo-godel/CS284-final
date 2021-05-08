#pragma once

#include "Snow.h"
#include "Particle.h"
#include "Vector2D.h"
#include "Constants.h"
#include "utilis.h"

#include <iostream>
#include <ctime>
#include <cstdlib>

using namespace std;

class Scene
{
public:
	vector<Snow*> snows;

	vector<Particle> particles;

	Scene() {};
	Scene(const Scene& s) : snows(s.snows) {}

	// init snow balls and convert them into particles
	void init();

	void draw();
};

