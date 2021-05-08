#include "Grid.h"
#include <cassert>

void Grid::initGridMassVel() {
	// Map particle to grid
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		// get the index of the grid cross point corresponding to the particle (it is on the bottom left of the particle)
		p.grid_p = (p.pos - origin) / node_size;
		assert(p.grid_p.x >= 0);
		assert(p.grid_p.y >= 0);
		int p_x = (int)p.grid_p.x;  // x coord index in grid
		int p_y = (int)p.grid_p.y;  // y coord index in grid

		// Map from (p_x - 1, p_y - 1) to (p_x + 2, p_y + 2)
		// The origin is bottom left, which means node_id = y * size.x + x
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			if (y < 0 || y >= size.y) // here size.y has already been added by 1
				continue;
			// Y interpolation
			double weight_y = bspline(p.grid_p.y - y);
			double dy = bsplineSlope(p.grid_p.y - y);

			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				if (x < 0 || x >= size.x)
					continue;

				// X interpolation
				double weight_x = bspline(p.grid_p.x - x);
				double dx = bsplineSlope(p.grid_p.x - x);

				// set weight of particles related nodes
				double w = weight_x * weight_y;
				p.weights[it] = w;

				// set weight gradient
				p.weight_gradient[it] = Vector2D(dx * weight_y, dy * weight_x);
				p.weight_gradient[it] /= node_size;

				// set node weighted mass and velocity
				int node_id = int(y * size.x + x);

				//nodes[node_id].mass += w * p.mass;
				atomicAdd(&(grid_ptr[node_id].mass), w * p.mass);
				//nodes[node_id].vel += p.vel * w * p.mass;
				atomicAdd(&(grid_ptr[node_id].vel), p.vel * w * p.mass);
				//nodes[node_id].active = true;
				atomicAdd(&(grid_ptr[node_id].active), 1);
			}
		}
	};

	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);

	thrust::for_each(
		thrust::device,
		nodes.begin(),
		nodes.end(),
		[=] __device__(Node& n) {
			if (n.active)
				n.vel /= n.mass;
		}
	);

}

// Calculate particles'volumes
void Grid::initVolumes() {
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		int p_x = (int)p.grid_p.x;
		int p_y = (int)p.grid_p.y;

		p.density = 0;
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				double w = p.weights[it];
				int node_id = int(y * size.x + x);
				if (w > BSPLINE_EPSILON) {
					p.density += w * grid_ptr[node_id].mass;
				}
			}
		}

		p.density /= node_area;
		p.volume = p.mass / p.density;
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
}

//Calculate grid's velocity of next timestep
void Grid::computeForce() {
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		// First calculate force based on mpmcourse
		Matrix2D I;
		Matrix2D R, S;
		Matrix2D U, Sig, V;

		polarDecomp(p.elastic_deformation, R, S);
		svd(p.elastic_deformation, U, Sig, V);

		double e = std::exp(HARDENING * (1.0f - p.plastic_deformation.det()));
		double lambda = LAMBDA * e;
		double mu = MU * e;
		double Je = Sig.det();

		Matrix2D temp = (p.elastic_deformation - U * V.T()) * p.elastic_deformation.T() * 2 * mu + I * lambda * Je * (Je - 1);
		temp = temp * p.volume;

		// accumulate particle stress to grids
		int p_x = (int)p.grid_p.x;
		int p_y = (int)p.grid_p.y;
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				double w = p.weights[it];
				int node_id = int(y * size.x + x);
				//Node& node = nodes[node_id];
				if (w > BSPLINE_EPSILON) {
					//grid_ptr[node_id].force -= temp * p.weight_gradient[it];
					atomicSub(&(grid_ptr[node_id].force), temp * p.weight_gradient[it]);
				}
			}
		}
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
}

void Grid::updateVelocity() {
    // here is how we update grid velocity
	thrust::for_each(
		thrust::device,
		nodes.begin(),
		nodes.end(),
		[=] __device__(Node& n) {
			if (n.active)
				n.vel_new = n.vel + TIMESTEP * (GRAVITY + n.force / n.mass);
		}
	);

    collisionGrid();
}

void Grid::updateDeformation() {
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		int p_x = (int)p.grid_p.x;
		int p_y = (int)p.grid_p.y;
		p.velocity_gradient = Matrix2D(0, 0, 0, 0);
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				double temp = p.weights[it];
				Vector2D delta_w = p.weight_gradient[it];
				int node_id = int(y * size.x + x);
				if (temp > BSPLINE_EPSILON) {
					p.velocity_gradient += outer_product(grid_ptr[node_id].vel_new, delta_w);
				}
			}
		}
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);

	auto func2 = [=] __device__(Particle & p) {
		Matrix2D I = Matrix2D();
		p.elastic_deformation = (I + p.velocity_gradient * TIMESTEP) * p.elastic_deformation;
		p.deformation_gradient = p.elastic_deformation * p.plastic_deformation;
		Matrix2D U, Sig, V;
		svd(p.elastic_deformation, U, Sig, V);
		for (int idx = 0; idx < 2; ++idx) {
			if (Sig(idx, idx) < CRIT_COMPRESS) {
				Sig(idx, idx) = CRIT_COMPRESS;
			}
			else if (Sig(idx, idx) > CRIT_STRETCH) {
				Sig(idx, idx) = CRIT_STRETCH;
			}
		}
		Matrix2D Sig_inv(1.0 / Sig(0, 0), 0, 0, 1.0 / Sig(1, 1));
		p.elastic_deformation = U * Sig * V.T();
		p.plastic_deformation = V * Sig_inv * U.T() * p.deformation_gradient;
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func2);
}

// Map back to particles
void Grid::updateParticlesVelocity() {
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		int p_x = (int)p.grid_p.x;
		int p_y = (int)p.grid_p.y;

		p.density = 0;

		Vector2D v_pic, v_flip = p.vel;
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				double w = p.weights[it];
				int node_id = int(y * size.x + x);
				if (w > BSPLINE_EPSILON) {
					//Node& node = nodes[node_id];
					v_pic += grid_ptr[node_id].vel_new * w;
					v_flip += (grid_ptr[node_id].vel_new - grid_ptr[node_id].vel) * w;
					p.density += w * grid_ptr[node_id].mass;
				}
			}
		}
		p.vel = v_flip * FLIP_PERCENT + v_pic * (1 - FLIP_PERCENT);
		p.density /= node_area;
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);

	collisionParticles();
}

void Grid::updateParticlesPosition() {
	auto func = [=] __device__(Particle & p) {
		p.pos += TIMESTEP * p.vel;
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
}

void Grid::collisionGrid() {

	auto func = [=] __device__(Node & n) {
		if (n.active) {
			Vector2D delta_scale = Vector2D(TIMESTEP, TIMESTEP);
			delta_scale /= node_size;
			Vector2D new_pos = n.vel_new * delta_scale + n.pos;
			if (new_pos.x < BSPLINE_RADIUS || new_pos.x > size.x - BSPLINE_RADIUS - 1) {
				n.vel_new.x = 0;
				n.vel_new.y *= STICKY;
			}
			if (new_pos.y < BSPLINE_RADIUS || new_pos.y > size.y - BSPLINE_RADIUS - 1) {
				n.vel_new.x *= STICKY;
				n.vel_new.y = 0;
			}
		}
	};
	thrust::for_each(thrust::device, nodes.begin(), nodes.end(), func);

}

void Grid::collisionParticles() {
	
	auto func = [=] __device__(Particle & p) {
		Vector2D delta_scale = Vector2D(TIMESTEP, TIMESTEP);
		delta_scale /= node_size;

		Vector2D new_pos = p.grid_p + p.vel * delta_scale;

		if (new_pos.x < BSPLINE_RADIUS - 1 || new_pos.x > size.x - BSPLINE_RADIUS) {
			p.vel.x *= -STICKY;
		}
		if (new_pos.y < BSPLINE_RADIUS - 1 || new_pos.y > size.y - BSPLINE_RADIUS) {
			p.vel.y *= -STICKY;
		}
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
}