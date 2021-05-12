#include "Scene.h"

double randomLR(double left, double right) {
	return left + rand() / (double)(RAND_MAX / (right - left));
}

void Scene::init() {
	// modify this part to generate anything interesting
	
	// First Scene
	/*
	Snow* snow1 = Snow::snowGenerator(Vector2D(1., .16), .15, Vector2D(0, 0));
	Snow* snow2 = Snow::snowGenerator(Vector2D(1., .4), .11, Vector2D(0, 0));
	Snow* snow3 = Snow::snowGenerator(Vector2D(0.5, .8), .1, Vector2D(15, -5));
	Snow* snow4 = Snow::snowGenerator(Vector2D(0.5, .61), .1, Vector2D(15, -5));
	Snow* snow5 = Snow::snowGenerator(Vector2D(1.5, .8), .1, Vector2D(-15, -5));
	Snow* snow6 = Snow::snowGenerator(Vector2D(1.5, .61), .1, Vector2D(-15, -5));
	*/
	

	// Second Scene
	
	Snow* snow1 = Snow::snowGenerator(Vector2D(1., .11), .1, Vector2D(0, 0));
	Snow* snow2 = Snow::snowGenerator(Vector2D(0.81, .11), .1, Vector2D(0, 0));
	Snow* snow3 = Snow::snowGenerator(Vector2D(1.19, .11), .1, Vector2D(0, 0));
	Snow* snow4 = Snow::snowGenerator(Vector2D(0.905, 0.2745), .1, Vector2D(0, 0));
	Snow* snow5 = Snow::snowGenerator(Vector2D(1.095, 0.2745), .1, Vector2D(0, 0));
	Snow* snow6 = Snow::snowGenerator(Vector2D(1, 0.439), .1, Vector2D(0, 0));
	Snow* snow7 = Snow::snowGenerator(Vector2D(1, 0.80), .15, Vector2D(0, -5));
	

	snows.push_back(snow1);
	snows.push_back(snow2);
	snows.push_back(snow3);
	snows.push_back(snow4);
	snows.push_back(snow5);
	snows.push_back(snow6);
	snows.push_back(snow7);


	// particles generator based on distribution of each shapes' area
	double total_area = 0;
	int len = snows.size();
	for (int i = 0; i < len; ++i) {
		total_area += snows[i]->area();
	}

	double p_area = PARTICLE_DIAM * PARTICLE_DIAM, p_mass = p_area * DENSITY;

	int total_num = total_area / p_area;
	particles.reserve(total_num);
	
	for (int i = 0, num_counts = 0; i < len; ++i) {
		Snow* s = snows[i];
		
		// propotional to area
		int temp_num = 0;
		if (i != len - 1)
			temp_num = s->area() / total_area * total_num;
		// final snow gets all remains
		else
			temp_num = total_num - num_counts;
		num_counts += temp_num;

		s->boundingBox();

		// randomly scatter particles into snow
		srand(time(0));

		while (temp_num) {
			double rx = randomLR(s->bounds[0], s->bounds[1]);
			double ry = randomLR(s->bounds[2], s->bounds[3]);

			//cout << rx << " " << ry << endl;

			if (s->ifInside(rx, ry)) {
				temp_num--;
				particles.push_back(Particle(Vector2D(rx, ry), Vector2D(s->vel), p_mass));
			}
		}
	}
}


void Scene::draw() {
	/*
	if (SUPPORTS_POINT_SMOOTH)
		glEnable(GL_POINT_SMOOTH);

	glPointSize(2);
	glBegin(GL_POINTS);

	int len = particles.size();
	for (int i = 0; i < len; i++) {
		Particle& p = particles[i];
		// We can use the particle's density to vary color
		float contrast = 0.6;
		float density = p.density / DENSITY * contrast;
		density += 1 - contrast;
		glColor3f(density, density, density);
		glVertex2f(p.pos.x, p.pos.y);
	}
	glEnd();
	*/
}