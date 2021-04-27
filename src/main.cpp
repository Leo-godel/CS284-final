#include <iostream>
#include <GLFW/glfw3.h>
#include "Constants.h"
#include "Scene.h"
#include "Simulator.h"

using namespace std;

Simulator* simluator;
Scene* scene;

void initGLContext() {
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, WIN_SIZE_X, WIN_SIZE_Y);
	glOrtho(0, WIN_METERS_X, 0, WIN_METERS_Y, 0, 1);
}

static void error_callback(int error, const char* description) {
	printf("\nError: %s", description);
}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	
}
void mouse_callback(GLFWwindow* window, int btn, int action, int mods) {
	
}

GLFWwindow* initGLFWContext() {
	glfwSetErrorCallback(error_callback);

	if (!glfwInit())
		exit(EXIT_FAILURE);

	GLFWwindow* window = glfwCreateWindow(WIN_SIZE_X, WIN_SIZE_Y, "Snow Simulator", NULL, NULL);

	if (!window) {
		glfwTerminate();
		exit(EXIT_FAILURE);
	}

	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_callback);

	return window;
}

void init() {
	scene = new Scene();
	scene->init();

	simluator = new Simulator(scene);
}

void update() {
	simluator->update();
}

void render() {
	scene->draw();
}

int main()
{
	GLFWwindow* window = initGLFWContext();
	initGLContext();

	
	init();

	int frame_count = 0;
	while (!glfwWindowShouldClose(window)) {
		cout << "frame:" << frame_count << endl;
		 update();

		glClearColor(0, 0, 0, 1);
		glClear(GL_COLOR_BUFFER_BIT);
		render();

		glfwSwapBuffers(window);
		glfwPollEvents();

		frame_count++;
	}

	//Exit
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
    return 0;
}