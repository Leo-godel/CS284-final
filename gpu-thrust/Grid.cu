#include "Grid.h"
#include <cassert>
#include <cstdio>

#define MATRIX_EPSILON 1e-6

__host__ __device__ double bspline(double x) {
	x = fabs(x);
	double w;
	if (x < 1)
		w = x * x * (x / 2 - 1) + 2 / 3.0;
	else if (x < 2)
		w = x * (x * (-x / 6 + 1) - 2) + 4 / 3.0;
	else return 0;

	return w;
}
//Slope of interpolation function
__host__ __device__ double bsplineSlope(double x) {
	double abs_x = fabs(x);
	if (abs_x < 1)
		return 1.5 * x * abs_x - 2 * x;
	else if (x < 2)
		return -x * abs_x / 2 + 2 * x - 2 * x / abs_x;
	else return 0;
}

__host__ __device__ Matrix2D outer_product(Vector2D& a, Vector2D& b) {
	return Matrix2D(a.x * b.x, a.x * b.y, a.y * b.x, a.y * b.y);
}

__host__ __device__ void polarDecomp(Matrix2D m, Matrix2D& R, Matrix2D& S) {
	auto x = m(0, 0) + m(1, 1);
	auto y = m(1, 0) - m(0, 1);
	auto scale = 1.0 / sqrt(x * x + y * y);
	auto c = x * scale, s = y * scale;
	R(0, 0) = c;
	R(0, 1) = -s;
	R(1, 0) = s;
	R(1, 1) = c;
	S = R.T() * m;
}

__host__ __device__ void svd(Matrix2D m, Matrix2D& U, Matrix2D& sig, Matrix2D& V){
    if (fabs(m(1, 0) - m(0, 1)) < MATRIX_EPSILON && fabs(m(1, 0)) < MATRIX_EPSILON) {
        U = Matrix2D(m(0, 0) < 0 ? -1 : 1, 0, 0, m(1, 1) < 0 ? -1 : 1);
        sig(0, 0) = fabs(m(0, 0)), sig(1, 1) = fabs(m(1, 1));
        V = Matrix2D();
    }
    else {
        double j = m(0, 0) * m(0, 0) + m(1, 0) * m(1, 0);
        double k = m(0, 1) * m(0, 1) + m(1, 1) * m(1, 1);
        double v_c = m(0, 0) * m(0, 1) + m(1, 0) * m(1, 1);

        if (fabs(v_c) < MATRIX_EPSILON) {
            double s1 = sqrt(j);
            double s2 = fabs(j - k) < MATRIX_EPSILON ? s1 : sqrt(k);
            sig(0, 0) = s1, sig(1, 1) = s2;
            V = Matrix2D();
            U = Matrix2D(m(0, 0) / s1, m(0, 1) / s2, m(1, 0) / s1, m(1, 1) / s2);
        }
        else {
            double jmk = j - k,
                    jpk = j + k,
                    root = sqrt(jmk * jmk + 4 * v_c * v_c),
                    eig = (jpk + root) / 2,
                    s1 = sqrt(eig),
                    s2 = fabs(root) < MATRIX_EPSILON ? s1 : sqrt((jpk - root) / 2);
            sig(0, 0) = s1, sig(1, 1) = s2;
            double v_s = eig - j,
                    len = sqrt(v_s * v_s + v_c * v_c);
            v_c /= len;
            v_s /= len;
            V = Matrix2D(v_c, -v_s, v_s, v_c);
            U = Matrix2D(
                    (m(0, 0) * v_c + m(0, 1) * v_s) / s1,
                    (m(0, 1) * v_c - m(0, 0) * v_s) / s2,
                    (m(1, 0) * v_c + m(1, 1) * v_s) / s1,
                    (m(1, 1) * v_c - m(1, 0) * v_s) / s2
            );
        }
    }
}

__device__ double my_atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                              (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}

struct OBCmp {
	__host__ __device__
	bool operator()(const SpGrid& o1, const SpGrid& o2) {
		return o1.node_id< o2.node_id;
	}
};

__host__ void Grid::initGridMassVel() {
	// Map particle to grid
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__ (Particle & p) {
		Vector2D _origin(0, 0), _node_size(1. / 128., 1. /128.), _size(256 + 1, 128 + 1);
		// get the index of the grid cross point corresponding to the particle (it is on the bottom left of the particle)
		p.grid_p = (p.pos - _origin) / _node_size;
		int p_x = (int)p.grid_p.x;  // x coord index in grid
		int p_y = (int)p.grid_p.y;  // y coord index in grid

		// Map from (p_x - 1, p_y - 1) to (p_x + 2, p_y + 2)
		// The origin is bottom left, which means node_id = y * size.x + x
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			if (y < 0 || y >= _size.y) // here size.y has already been added by 1
				continue;
			// Y interpolation
			double weight_y = bspline(p.grid_p.y - y);
			double dy = bsplineSlope(p.grid_p.y - y);

			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				if (x < 0 || x >= _size.x)
					continue;

				// X interpolation
				double weight_x = bspline(p.grid_p.x - x);
				double dx = bsplineSlope(p.grid_p.x - x);

				// set weight of particles related nodes
				double w = weight_x * weight_y;
				p.weights[it] = w;

				// set weight gradient
				p.weight_gradient[it] = Vector2D(dx * weight_y, dy * weight_x);
				p.weight_gradient[it] /= _node_size;

				// set node weighted mass and velocity
				int node_id = int(y * _size.x + x);

				//nodes[node_id].mass += w * p.mass;
				//my_atomicAdd(&(grid_ptr[node_id].mass), w * p.mass);
				//nodes[node_id].vel += p.vel * w * p.mass;
				Vector2D temp = p.vel * w * p.mass;
				//my_atomicAdd(&(grid_ptr[node_id].vel.x), temp.x);
				//my_atomicAdd(&(grid_ptr[node_id].vel.y), temp.y);
				//nodes[node_id].active = true;
				//atomicAdd(&(grid_ptr[node_id].active), 1);
			}
		}
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
	//cout << other_particles.size() << endl;
	thrust::sort(thrust::device, other_particles.begin(), other_particles.end(), OBCmp());

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
__host__ void Grid::initVolumes() {
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		Vector2D _origin(0, 0), _node_size(1. / 128., 1. /128.), _size(256 + 1, 128 + 1);
		double _node_area = _node_size.x * _node_size.y;

		int p_x = (int)p.grid_p.x;
		int p_y = (int)p.grid_p.y;

		p.density = 0;
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				if (y < 0 || y >= _size.y || x < 0 || x >= _size.x)
					continue;
				
				double w = p.weights[it];
				int node_id = int(y * _size.x + x);
				if (w > BSPLINE_EPSILON) {
					p.density += w * grid_ptr[node_id].mass;
				}
			}
		}

		p.density /= _node_area;
		p.volume = p.mass / p.density;
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
}

//Calculate grid's velocity of next timestep
__host__ void Grid::computeForce() {
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		// First calculate force based on mpmcourse
		Vector2D _origin(0, 0), _node_size(1. / 128., 1. /128.), _size(256 + 1, 128 + 1);
		double _node_area = _node_size.x * _node_size.y;

		Matrix2D I;
		Matrix2D U, Sig, V;

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
				if (y < 0 || y >= _size.y || x < 0 || x >= _size.x)
					continue;

				double w = p.weights[it];
				int node_id = int(y * _size.x + x);
				
				//Node& node = nodes[node_id];
				if (w > BSPLINE_EPSILON) {
					//grid_ptr[node_id].force -= temp * p.weight_gradient[it];
					Vector2D value = temp * p.weight_gradient[it];
					//my_atomicAdd(&(grid_ptr[node_id].force.x), -value.x);
					//my_atomicAdd(&(grid_ptr[node_id].force.y), -value.y);
				}
			}
		}
	};
	thrust::sort(thrust::device, other_particles.begin(), other_particles.end(), OBCmp());
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
}

__host__ void Grid::updateVelocity() {
    // here is how we update grid velocity
	thrust::for_each(
		thrust::device,
		nodes.begin(),
		nodes.end(),
		[=] __device__(Node& n) {
			double timestep = 0.0001;
			Vector2D gravity(0, -9.8);

			if (n.active) {
				n.vel_new = n.vel + timestep * (gravity + n.force / n.mass);
				//printf("updated!\n");
			}
		}
	);
    collisionGrid();
}

__host__ void Grid::updateDeformation() {
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		Vector2D _origin(0, 0), _node_size(1. / 128., 1. /128.), _size(256 + 1, 128 + 1);
		double _node_area = _node_size.x * _node_size.y;

		int p_x = (int)p.grid_p.x;
		int p_y = (int)p.grid_p.y;
		p.velocity_gradient = Matrix2D(0, 0, 0, 0);
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				if (y < 0 || y >= _size.y || x < 0 || x >= _size.x)
					continue;
				
				double temp = p.weights[it];
				Vector2D delta_w = p.weight_gradient[it];
				int node_id = int(y * _size.x + x);
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
__host__ void Grid::updateParticlesVelocity() {
	Node* grid_ptr = thrust::raw_pointer_cast(&nodes[0]);

	auto func = [=] __device__(Particle & p) {
		Vector2D _origin(0, 0), _node_size(1. / 128., 1. /128.), _size(256 + 1, 128 + 1);
		double _node_area = _node_size.x * _node_size.y;

		int p_x = (int)p.grid_p.x;
		int p_y = (int)p.grid_p.y;

		p.density = 0;

		Vector2D v_pic, v_flip = p.vel;
		for (int it = 0, y = p_y - 1; y <= p_y + 2; ++y) {
			for (int x = p_x - 1; x <= p_x + 2; ++x, ++it) {
				if (y < 0 || y >= _size.y || x < 0 || x >= _size.x)
					continue;
				
				double w = p.weights[it];
				int node_id = int(y * _size.x + x);
				if (w > BSPLINE_EPSILON) {
					//Node& node = nodes[node_id];
					v_pic += grid_ptr[node_id].vel_new * w;
					v_flip += (grid_ptr[node_id].vel_new - grid_ptr[node_id].vel) * w;
					p.density += w * grid_ptr[node_id].mass;
				}
			}
		}
		double flip_percent = .95;
		p.vel = v_flip * flip_percent + v_pic * (1 - flip_percent);
		p.density /= _node_area;
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);

	collisionParticles();
}

__host__ void Grid::updateParticlesPosition() {
	auto func = [=] __device__(Particle & p) {
		double timestep = 0.0001;
		p.pos += timestep * p.vel;
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
}

__host__ void Grid::collisionGrid() {

	auto func = [=] __device__(Node & n) {
		Vector2D _origin(0, 0), _node_size(1. / 128., 1. /128.), _size(256 + 1, 128 + 1);
		double _node_area = _node_size.x * _node_size.y;
		double timestep = 0.0001;

		if (n.active) {
			Vector2D delta_scale = Vector2D(timestep, timestep);
			delta_scale /= _node_size;
			Vector2D new_pos = n.vel_new * delta_scale + n.pos;
			if (new_pos.x < BSPLINE_RADIUS || new_pos.x > _size.x - BSPLINE_RADIUS - 1) {
				n.vel_new.x = 0;
				n.vel_new.y *= STICKY;
			}
			if (new_pos.y < BSPLINE_RADIUS || new_pos.y > _size.y - BSPLINE_RADIUS - 1) {
				n.vel_new.x *= STICKY;
				n.vel_new.y = 0;
			}
		}
	};
	thrust::for_each(thrust::device, nodes.begin(), nodes.end(), func);

}

__host__ void Grid::collisionParticles() {
	
	auto func = [=] __device__(Particle & p) {
		Vector2D _origin(0, 0), _node_size(1. / 128., 1. /128.), _size(256 + 1, 128 + 1);
		double _node_area = _node_size.x * _node_size.y;
		double timestep = 0.0001;

		Vector2D delta_scale = Vector2D(timestep, timestep);
		delta_scale /= _node_size;

		Vector2D new_pos = p.grid_p + p.vel * delta_scale;

		if (new_pos.x < BSPLINE_RADIUS - 1 || new_pos.x > _size.x - BSPLINE_RADIUS) {
			p.vel.x *= -STICKY;
		}
		if (new_pos.y < BSPLINE_RADIUS - 1 || new_pos.y > _size.y - BSPLINE_RADIUS) {
			p.vel.y *= -STICKY;
		}
	};
	thrust::for_each(thrust::device, particles.begin(), particles.end(), func);
}
