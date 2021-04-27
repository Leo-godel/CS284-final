#include "Grid.h"

void Grid::initGridMassVel() {
	//clear the nodes
	memset(nodes, 0, sizeof(Node) * nodes_length);

	// Map particle to grid
	int len = particles->size();
	for (int i = 0; i < len; ++i) {
		Particle& p = particles->at(i);

		// get the index of the grid cross point corresponding to the particle (it is on the bottom left of the particle)
		p.grid_p = (p.pos - origin) / node_size;
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
				// PS: Dont know why need this, reserve in case for particle inner force calculation
				p.weight_gradient[it] = Vector2D(dx * weight_y, dy * weight_x);
				p.weight_gradient[it] /= node_size;

				// set node weighted mass and velocity
				int node_id = int(y * size.x + x);
				nodes[node_id].mass += w * p.mass;
                nodes[node_id].vel += p.vel * w * p.mass;
			}
		}
	}

    for (int i = 0; i < nodes_length; ++i) {
        if(nodes[i].active)
            nodes[i].vel /= nodes[i].mass;
    }
}

// Map particles'velocity to grid
void Grid::initGridVel() {
	int len = particles->size();
	for (int i = 0; i < len; ++i) {
		Particle& p = particles->at(i);
		double p_x = p.grid_p.x, p_y = p.grid_p.y;
		
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				double w = p.weights[it];
				int node_id = int(y * size.x + x);
				if (w > BSPLINE_EPSILON) {
					nodes[node_id].vel += p.vel * w * p.mass;
					nodes[node_id].active = true;
				}
			}
		}
	}

	for (int i = 0; i < nodes_length; ++i) {
		if (nodes[i].active) {
			nodes[i].vel /= nodes[i].mass;
		}
	}

	collisionGrid();
}

// Calculate particles'volumes
void Grid::initVolumes() {
	int len = particles->size();

	for (int i = 0; i < len; ++i) {
		Particle& p = particles->at(i);
		int p_x = (int)p.grid_p.x;
		int p_y = (int)p.grid_p.y;

		p.density = 0;
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				double w = p.weights[it];
				int node_id = int(y * size.x + x);
				if (w > BSPLINE_EPSILON) {
					p.density += w * nodes[node_id].mass;
				}
			}
		}

		p.density /= node_area;
		p.volume = p.mass / p.density;
	}
}

//Calculate grid's velocity of next timestep
void Grid::computeForce() {
    int len = particles->size();

    for (int i = 0; i < len; ++i) {
        Particle &p = particles->at(i);

        // First calculate force based on mpmcourse and taichi 88 line
        Matrix2D R, S;
        polarDecomp(p.deformation_gradient, R, S);
        p.volume = p.deformation_gradient.det();

        double e = std::exp(HARDENING * (1.0f - p.volume));
        auto lambda = LAMBDA * e;
        auto mu = MU * e;

        Matrix2D F_inv = p.deformation_gradient.inv();
        Matrix2D P = ((p.deformation_gradient - R) * mu * 2) + (F_inv.T() * lambda * (p.volume - 1) * p.volume);
        p.stress = P * p.deformation_gradient.T() * (1 / p.volume);

        // accumulate particle stress to grids
        int p_x = (int)p.grid_p.x;
        int p_y = (int)p.grid_p.y;
        for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
            for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
                double w = p.weights[it];
                int node_id = int(y * size.x + x);
                if (w > BSPLINE_EPSILON) {
                    nodes[node_id].force -= p.stress * p.weight_gradient[it] * p.volume;
                }
            }
        }
    }
    collisionGrid();
}

void Grid::updateVelocity() {
    // here is how we update grid velocity
    for (int i = 0; i < nodes_length; ++i) {
        if (nodes[i].active) {
            nodes[i].vel_new = nodes[i].vel + TIMESTEP * (GRAVITY + nodes[i].force / nodes[i].mass);
        }
    }
}

void Grid::updateDeformation() {
    int len = particles->size();

    for (int i = 0; i < len; ++i) {
        Particle& p = particles->at(i);
        int p_x = (int)p.grid_p.x;
        int p_y = (int)p.grid_p.y;
        for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
            for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
                double temp = p.weights[it];
                Vector2D delta_w = p.weight_gradient[it];
                int node_id = int(y * size.x + x);
                if (temp > BSPLINE_EPSILON) {
                    p.velocity_gradient += outer_product(nodes[node_id].vel, delta_w);
                }
            }
        }
    }
    Matrix2D I = Matrix2D();
    for(int i = 0; i < len; ++i) {
        Particle& p = particles->at(i);
        p.deformation_gradient = (I + p.velocity_gradient * TIMESTEP) * p.deformation_gradient;
    }
}

// Map back to particles
void Grid::updateParticlesVelocity() {
	int len = particles->size();

	for (int i = 0; i < len; ++i) {
		Particle& p = particles->at(i);
        int p_x = (int)p.grid_p.x;
        int p_y = (int)p.grid_p.y;

		Vector2D v_pic, v_flip = p.vel;

		p.density = 0;
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				double w = p.weights[it];
				int node_id = int(y * size.x + x);
				if (w > BSPLINE_EPSILON) {
					Node& node = nodes[node_id];
					v_pic += node.vel_new * w;
					v_flip += (node.vel_new - node.vel) * w;
				}
			}
		}
		p.vel = v_flip * FLIP_PERCENT + v_pic * (1 - FLIP_PERCENT);
	}
	collisionParticles();
}

void Grid::collisionGrid() {
	Vector2D delta_scale = Vector2D(TIMESTEP, TIMESTEP);
	delta_scale /= node_size;

	for (int it = 0, y = 0; y < (int)size.y; ++y) {
		for (int x = 0; x < (int)size.x; ++x, ++it) {
			Node& node = nodes[it];
			if (!node.active)
				continue;

			Vector2D new_pos = node.vel_new * delta_scale + Vector2D(x, y);
			if (new_pos.x < BSPLINE_RADIUS || new_pos.x > size.x - BSPLINE_RADIUS - 1) {
				node.vel_new.x = 0;
				node.vel_new.y *= STICKY;
			}
			if (new_pos.y < BSPLINE_RADIUS || new_pos.y > size.y - BSPLINE_RADIUS - 1) {
				node.vel_new.x *= STICKY;
				node.vel_new.y = 0;
			}
		}
	}
}

void Grid::collisionParticles() {
	Vector2D delta_scale = Vector2D(TIMESTEP, TIMESTEP);
	delta_scale /= node_size;
	int len = particles->size();

	for (int i = 0; i < len; ++i) {
		Particle& p = particles->at(i);
		
		Vector2D new_pos = p.grid_p + p.vel * delta_scale;

		if (new_pos.x < BSPLINE_RADIUS-1 || new_pos.x > size.x - BSPLINE_RADIUS) {
			p.vel.x *= -STICKY;
		}
		if (new_pos.y < BSPLINE_RADIUS-1 || new_pos.y > size.y - BSPLINE_RADIUS) {
			p.vel.y *= -STICKY;
		}
	}
}