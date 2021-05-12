#ifndef GRID_H
#define GRID_H

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
#include <thrust/copy.h>
#include <thrust/device_vector.h>

struct Node {
	double mass;
	int active;
	Vector2D pos, vel, vel_new, force;
	__host__ __device__ Node() {
		mass = 0;
		active = 0;
		pos = vel = vel_new = force = Vector2D();
	}
};

struct SpGrid {
	int node_id, part_id;
	__host__ __device__ SpGrid() {
		node_id = part_id = 0;
	}
};

class Grid
{
public:
	// Grid origin and size
	Vector2D origin, size;
	
	// Grid nodes
	Vector2D node_size;
	double node_area;
	int nodes_length;

	__host__ Grid(Vector2D pos, Vector2D dims, Vector2D window_size, Scene* s) :
		origin(pos), size(window_size) {
		//particles = &s->particles;

		node_size = dims / window_size;
		size.x += 1, size.y += 1;

		nodes_length = int(size.x * size.y);
		//nodes = new Node[nodes_length];
		node_area = node_size.x * node_size.y;

		cout << "# particles: " << s->particles.size() << endl;

		particles.resize(s->particles.size());
		thrust::copy(s->particles.begin(), s->particles.end(), particles.begin());

		
		other_particles.resize(16 * s->particles.size());
		/*
		for (int i = 0; i < 16; ++i) {
			int len = i * s->particles.size();
			thrust::copy(s->particles.begin(), s->particles.end(), other_particles.begin() + len);
		}*/

		nodes.resize(nodes_length);
		vector<Node> temp_nodes;
		temp_nodes.resize(nodes_length);
		for (int y = 0; y < size.y; ++y) {
			for (int x = 0; x < size.x; ++x) {
				int node_id = int(y * size.x + x);
				temp_nodes[node_id].pos = Vector2D(x, y);
			}
		}
		thrust::copy(temp_nodes.begin(), temp_nodes.end(), nodes.begin());
	}

	// Map particles to Grid
	__host__ void initGridMassVel();
	// Map volumes of particles: contains how many nodes
	// only called once
	__host__ void initVolumes();
	// Compute force and velocity of next time step
	__host__ void computeForce();
    __host__ void updateVelocity();
    __host__ void updateDeformation();
	// Map back to pariticles
	__host__ void updateParticlesVelocity();
	__host__ void updateParticlesPosition();

	// Board Collision
	__host__ void collisionGrid();
	__host__ void collisionParticles();

private:
	thrust::device_vector<Particle> particles;
	thrust::device_vector<SpGrid> other_particles;
	thrust::device_vector<Node> nodes;
};

#endif