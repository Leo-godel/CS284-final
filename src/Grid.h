#pragma once

#include "Constants.h"
#include "Vector2D.h"
#include "Scene.h"

#include <cmath>
#include <cstring>

struct Node {
	double mass;
	bool active;
	Vector2D vel, vel_new;
};

class Grid
{
public:
	vector<Particle>* particles;
	// Grid origin and size
	Vector2D origin, size;
	
	// Grid nodes
	Node* nodes;
	Vector2D node_size;
	double node_area;
	int nodes_length;

	Grid(Vector2D pos, Vector2D dims, Vector2D window_size, Scene* s) :
		origin(pos), size(window_size) {
		particles = &s->particles;

		node_size = dims / window_size;
		size.x += 1, size.y += 1;

		nodes_length = int(size.x * size.y);
		nodes = new Node[nodes_length];
		node_area = node_size.x * node_size.y;
	}

	// Map particles to Grid
	void initGridMassVel();
	void initGridVel();
	// Map volumes of particles: contains how many nodes
	// only called once
	void initVolumes();
	// Compute velocity of next time step
	void computeForceVel();
	// Map back to pariticles
	void updateParticles();

	// Board Collision
	void collisionGrid();
	void collisionParticles();

	void draw();
};

