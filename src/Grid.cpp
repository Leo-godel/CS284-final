#include "Grid.h"
#include <cassert>

void Grid::initGridMassVel() {
	//clear the nodes
	memset(nodes, 0, sizeof(Node) * nodes_length);

	// Map particle to grid
	int len = particles->size();
	for (int i = 0; i < len; ++i) {
		Particle& p = particles->at(i);

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
				if (!nodes[node_id].active) {
				    assert(nodes[node_id].mass == 0);
                    assert(nodes[node_id].vel.x == 0);
                    assert(nodes[node_id].vel.y == 0);
				}
				nodes[node_id].mass += w * p.mass;
                nodes[node_id].vel += p.vel * w * p.mass;
                nodes[node_id].active = true;
			}
		}
	}

    for (int i = 0; i < nodes_length; ++i) {
        if(nodes[i].active)
            nodes[i].vel /= nodes[i].mass;
        else {
            assert(nodes[i].mass == 0);
            assert(nodes[i].vel.x == 0);
            assert(nodes[i].vel.y == 0);
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

        // First calculate force based on mpmcourse
        Matrix2D I;
        Matrix2D R, S;
        Matrix2D U, Sig, V;
//        assert(p.deformation_gradient(0,0) - (p.elastic_deformation * p.plastic_deformation)(0, 0) < 1e-6);
//        assert(p.deformation_gradient(1,0) - (p.elastic_deformation * p.plastic_deformation)(1, 0) < 1e-6);
//        assert(p.deformation_gradient.det() - (p.elastic_deformation * p.plastic_deformation).det()  < 1e-6 );
        polarDecomp(p.elastic_deformation, R, S);
        svd(p.elastic_deformation, U, Sig, V);
//        assert(p.elastic_deformation.det() -  (U * Sig * V.T()).det() < 1e-6);

        double e = std::exp(HARDENING * (1.0f - p.plastic_deformation.det()));
        double lambda = LAMBDA * e;
        double mu = MU * e;
        double Je = Sig.det();
//        std::cout << "Je: " << Je << endl;
//        std::cout << "harden: " << e << std::endl;

//        Matrix2D F_inv = p.deformation_gradient.inv();
//        Matrix2D P = ((p.deformation_gradient - R) * mu * 2) + (F_inv.T() * lambda * (p.volume - 1) * p.volume);
//        p.stress = P * p.deformation_gradient.T() * (1 / p.deformation_gradient.det());
//        cout << "stress: " << p.stress.det() << endl;
//        assert(p.deformation_gradient(0,0) == p.elastic_deformation(0, 0));
//        std::cout << "elastic: " << p.elastic_deformation(0, 0) << " " << p.elastic_deformation(0, 1) << " " <<
//        p.elastic_deformation(1, 0) << " " << p.elastic_deformation(1, 1) << std::endl;
//        double Je = Sig(0,0) * Sig(1,1);
//        assert(Je == p.deformation_gradient.det());
//        Matrix2D temp = (p.elastic_deformation - U * V.T()) * p.elastic_deformation.T() * 2 * mu + I * lambda * (Je - 1) * Je;
//        std::cout << "lambda: " << lambda << "mu: " << mu << endl;
//        std::cout << "Je: " << Je << endl;
//        std::cout << "temp: " << temp(0, 0) << " " << temp(0, 1) << " " << temp(1, 0) << " " << temp(1, 1) << std::endl;
//        std::cout << "R: " << R(0, 0) << " " << R(0, 1) << " " << R(1, 0) << " " << R(1, 1) << std::endl;
//        p.stress = temp * (1 / p.deformation_gradient.det());

//        cout << "comp1: " << ((p.elastic_deformation - R) * p.elastic_deformation.T()).det() << endl;
//        cout << "elastic: " << p.elastic_deformation.det() << endl;
//        cout << "stress: " << p.stress.det() << endl;

        Matrix2D temp = (p.elastic_deformation - U * V.T()) * p.elastic_deformation.T() * 2 * mu + I * lambda * Je * (Je - 1);
		temp = temp * p.volume;
//        std::cout << "temp: " << temp(0, 0) << " " << temp(1, 0) << " " << temp(0, 1) << " " << temp(1, 1) << endl;
//        std::cout << "temp det: " << temp.det() << endl;

        // accumulate particle stress to grids
        int p_x = (int)p.grid_p.x;
        int p_y = (int)p.grid_p.y;
        for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
            for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
                double w = p.weights[it];
                int node_id = int(y * size.x + x);
                Node& node = nodes[node_id];
                if (w > BSPLINE_EPSILON) {
//                    std::cout << p.volume << endl;
                    node.force -= temp * p.weight_gradient[it];
                }
            }
        }
    }
}

void Grid::updateVelocity() {
    // here is how we update grid velocity
    for (int i = 0; i < nodes_length; ++i) {
        if (nodes[i].active) {
//            cout << "force: " << (nodes[i].force).norm() << endl;
//            cout << "acceration: " << (nodes[i].force / nodes[i].mass).norm() << endl;
            nodes[i].vel_new = nodes[i].vel + TIMESTEP * (GRAVITY + nodes[i].force / nodes[i].mass);
        }
    }
    collisionGrid();
}

void Grid::updateDeformation() {
    int len = particles->size();

    for (int i = 0; i < len; ++i) {
        Particle& p = particles->at(i);
        int p_x = (int)p.grid_p.x;
        int p_y = (int)p.grid_p.y;
        p.velocity_gradient = Matrix2D(0, 0, 0, 0);
        for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
            for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
                double temp = p.weights[it];
                Vector2D delta_w = p.weight_gradient[it];
                int node_id = int(y * size.x + x);
                if (temp > BSPLINE_EPSILON) {
//                    std::cout << "Weight gradient: " << delta_w[0] << " " << delta_w[1] << std::endl;
                    p.velocity_gradient += outer_product(nodes[node_id].vel_new, delta_w);
                }
            }
        }
    }

    Matrix2D I = Matrix2D();
    for(int i = 0; i < len; ++i) {
        Particle& p = particles->at(i);
//        p.deformation_gradient = (I + p.velocity_gradient * TIMESTEP) * p.deformation_gradient;
//        std::cout << p.velocity_gradient(0, 0) << " " << p.velocity_gradient(0, 1) << endl;
//        std::cout << "Velocity gradient: " << p.velocity_gradient.det() << endl;
//        std::cout << "deformation update ratio: " << (I + p.velocity_gradient * TIMESTEP).det() << endl;
        p.elastic_deformation = (I + p.velocity_gradient * TIMESTEP) * p.elastic_deformation;
        p.deformation_gradient = p.elastic_deformation * p.plastic_deformation;
        Matrix2D U, Sig, V;
        svd(p.elastic_deformation, U, Sig, V);
        for (int idx = 0; idx < 2; ++idx){
            if (Sig(idx, idx) < CRIT_COMPRESS) {
                Sig(idx, idx) = CRIT_COMPRESS;
//                cout << "Clip compress. " << endl;
            }
            else if (Sig(idx, idx) > CRIT_STRETCH) {
//                cout << "Clip stretch. " << endl;
                Sig(idx, idx) = CRIT_STRETCH;
            }
        }
        Matrix2D Sig_inv(1.0 / Sig(0, 0), 0, 0, 1.0 / Sig(1, 1));
        p.elastic_deformation = U * Sig * V.T();
        p.plastic_deformation = V * Sig_inv * U.T() * p.deformation_gradient;
    }
}

// Map back to particles
void Grid::updateParticlesVelocity() {
	int len = particles->size();

	for (int i = 0; i < len; ++i) {
		Particle& p = particles->at(i);
        int p_x = (int)p.grid_p.x;
        int p_y = (int)p.grid_p.y;

        p.density = 0;

		Vector2D v_pic, v_flip = p.vel;
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				double w = p.weights[it];
				int node_id = int(y * size.x + x);
				if (w > BSPLINE_EPSILON) {
					Node& node = nodes[node_id];
					v_pic += node.vel_new * w;
					v_flip += (node.vel_new - node.vel) * w;
                    p.density += w * node.mass;
				}
			}
		}
		p.vel = v_flip * FLIP_PERCENT + v_pic * (1 - FLIP_PERCENT);
        p.density /= node_area;
	}
	collisionParticles();
}

void Grid::updateParticlesPosition() {
	int len = particles->size();
	for (int i = 0; i < len; ++i) {
		Particle& p = particles->at(i);
		
		p.pos += TIMESTEP * p.vel;
	}
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