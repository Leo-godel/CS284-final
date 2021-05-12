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
		node_area = node_size.x * node_size.y;

		particle_count = s->particles.size();
		cout << "parts size: " << particle_count << endl;
		cout << "part per thread: " << particle_count / 128 / 128 << endl;
		particles = new Particle[particle_count];
		nodes = new Node[nodes_length];

		for(int i = 0; i < particle_count; ++i) {
			particles[i] = s->particles[i];
		}

    	cudaMalloc((void**)&particles_gpu, particle_count * sizeof(Particle));
    	cudaMemcpy(particles_gpu, particles, particle_count * sizeof(Particle), cudaMemcpyHostToDevice);
		
		for (int y = 0; y < size.y; ++y) {
			for (int x = 0; x < size.x; ++x) {
				int node_id = int(y * size.x + x);
				nodes[node_id].pos = Vector2D(x, y);
			}
		}

		cudaMalloc((void**)&nodes_gpu, nodes_length * sizeof(Node));
    	cudaMemcpy(nodes_gpu, nodes, nodes_length * sizeof(Node), cudaMemcpyHostToDevice);
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
	__host__ void G2P2G();
	// Map back to pariticles
	__host__ void updateParticlesVelocity();
	__host__ void updateParticlesPosition();

	// Board Collision
	__host__ void collisionGrid();
	__host__ void collisionParticles();

private:
	// thrust::device_vector<Particle> particles;
	// thrust::device_vector<Node> nodes;
	Particle* particles;
	Particle* particles_gpu;
	int particle_count;
	Node* nodes;
	Node* nodes_gpu;
};

#endif