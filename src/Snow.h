#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <GLFW/glfw3.h>
#include "Vector2D.h"

using namespace std;

class Snow
{
public:
	vector<Vector2D> vertices;
	Vector2D vel;
	double bounds[4];

	Snow(Vector2D vel = Vector2D()) : vel(vel) {}
	Snow(const Snow& s) : vertices(s.vertices), vel(s.vel) {}

	//add vertices
	void addVertices(double x, double y);

	//check inside
	bool ifInside(double x, double y);

	// area
	double area();

	// bounding box
	void boundingBox();

	// ball shape snow generator
	static Snow* snowGenerator(Vector2D o, double radius, Vector2D vel);

};

