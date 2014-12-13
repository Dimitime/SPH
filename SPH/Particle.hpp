/*
 * Class definition for a Particle
 *
 */

#include <glm/glm.hpp>
#include <GL/freeglut.h>

class Particle {
public:
    //Default constructor
	Particle() {
        pos = glm::vec3(0.0f,0.0f,0.0f);
        vel = glm::vec3(0.0f,0.0f,0.0f);
        accel = glm::vec3(0.0f,0.0f,0.0f);

        mass = 0.0f;
        density = 0.0f;
        pressure = 0.0f;
        thermal = 0.0f;
    }

    //Explicit constructor
	Particle(glm::vec3 npos, glm::vec3 nvel, glm::vec3 naccel, float nmass, float ndensity, float npressure, float nthermal) {
        pos = npos;
        vel = nvel;
        accel = naccel;
        mass = nmass;
        density = ndensity;
        pressure = npressure;
        thermal = nthermal;
    }

    //Destructor
	~Particle(){}

    //Displays this particle
    //TODO use a real sphere making algorithm, not this crap
	void display(){
        glutSolidSphere(1, 10, 10);
    }

//private:
    glm::vec3 pos;
    glm::vec3 vel;
    glm::vec3 accel;

    float mass;
    float density;
    float pressure;
    float thermal;
};
