/*
 * Class definition for a Sphere fairly straightforward.
 */

#include <glm/glm.hpp>

class Sphere {
public:
    //Default constructor
	Sphere() {
        pos = glm::vec3(0.0f,0.0f,0.0f);
        vel = glm::vec3(0.0f,0.0f,0.0f);
        radius = 1.0f;
		mass = 1.0f;
    }

    //Explicit constructor
	Sphere(glm::vec3 npos, glm::vec3 nvel, float nradius, float nmass) {
		pos = npos;
		vel = nvel;
		radius = nradius;
		mass = nmass;
    }

    //Destructor
	~Sphere(){}

//private:
    glm::vec3 pos;
    glm::vec3 vel;
    glm::vec3 accel;
	float radius;
	float mass;
};