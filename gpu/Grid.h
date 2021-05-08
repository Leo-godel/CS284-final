#pragma once

#include "Constants.h"
#include "Vector2D.h"
#include "Scene.h"

#include <cmath>
#include <cstring>
#include <vector>

#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include <thrust/copy.h>
#include <thrust/device_vector.h>


struct Node {
	double mass;
	int active;
	Vector2D pos, vel, vel_new, force;
	Node() {
		mass = 0;
		active = 0;
		pos = vel = vel_new = force = Vector2D();
	}
};

class Grid
{
public:
	thrust::device_vector<Particle> particles;
	thrust::device_vector<Node> nodes;
	// Grid origin and size
	Vector2D origin, size;
	
	// Grid nodes
	Vector2D node_size;
	double node_area;
	int nodes_length;

	Grid(Vector2D pos, Vector2D dims, Vector2D window_size, Scene* s) :
		origin(pos), size(window_size) {
		//particles = &s->particles;

		node_size = dims / window_size;
		size.x += 1, size.y += 1;

		nodes_length = int(size.x * size.y);
		//nodes = new Node[nodes_length];
		node_area = node_size.x * node_size.y;

		particles.resize(s->particles.size());
		nodes.resize(nodes_length);
		Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);
		for (int y = 0; y < size.y; ++y) {
			for (int x = 0; x < size.x; ++x) {
				int node_id = int(y * size.x + x);
				grid_ptr[node_id].pos = Vector2D(x, y);
			}
		}
	}

	// Map particles to Grid
	void initGridMassVel();
	// Map volumes of particles: contains how many nodes
	// only called once
	void initVolumes();
	// Compute force and velocity of next time step
	void computeForce();
    void updateVelocity();
    void updateDeformation();
	// Map back to pariticles
	void updateParticlesVelocity();
	void updateParticlesPosition();

	// Board Collision
	void collisionGrid();
	void collisionParticles();

	void draw();
};

