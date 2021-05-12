#include "Snow.h"

void Snow::addVertices(double x, double y) {
	vertices.push_back(Vector2D(x, y));
}

//Source: http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html	
bool Snow::ifInside(double x, double y) {
	bool res = false;
	int len = vertices.size();
	for (int i = 0, j = len - 1; i < len; j = i++) {
		Vector2D &vi = vertices[i], &vj = vertices[j];
		bool condition1 = (vi[1] > y) != (vj[1] > y);
		bool condition2 = x < (vj[0] - vi[0]) * (y - vi[1]) / (vj[1] - vi[1]) + vi[0];
		if (condition1 && condition2)
			res = !res;
	}
	return res;
}

//Source: http://www.mathopenref.com/coordpolygonarea2.html
double Snow::area() {
	double area = 0;
	int len = vertices.size();
	for (int i = 0, j = len - 1; i < len; j = i++) {
		Vector2D &vi = vertices[i], &vj = vertices[j];
		area += (vj[0] + vi[0]) * (vj[1] - vi[1]);
	}
	return fabs(area / 2);
}

void Snow::boundingBox() {
	// bounds[0] = min_x, bounds[1] = max_x
	bounds[0] = bounds[1] = vertices[0].x;
	// bounds[2] = min_y, bounds[3] = max_y
	bounds[2] = bounds[3] = vertices[0].y;

	int len = vertices.size();
	for (int i = 0; i < len; ++i) {
		Vector2D& vi = vertices[i];

		bounds[0] = fmin(bounds[0], vi.x);
		bounds[1] = fmax(bounds[1], vi.x);
		bounds[2] = fmin(bounds[2], vi.y);
		bounds[3] = fmax(bounds[3], vi.y);
	}
}

//Cool circle algorithm: http://slabode.exofire.net/circle_draw.shtml
Snow* Snow::snowGenerator(Vector2D o, double radius, Vector2D vel) {
	Snow* snow = new Snow(vel);
	
	// sides number of snow
	const int sides = 18;
	const double theta = 6.283185307 / (double)sides;
	const double tan_fac = tan(theta), cos_fac = cos(theta);

	double x = radius, y = 0;

	for (int i = 0; i < sides; ++i) {
		snow->addVertices(x + o.x, y + o.y);
		double flip_x = -y, flip_y = x;
		x += flip_x * tan_fac;
		y += flip_y * tan_fac;
		x *= cos_fac;
		y *= cos_fac;
	}
	return snow;
}