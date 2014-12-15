/*
 * Class definition for a Particle
 *
 */

#include "Particle.hpp"
#include "Wall.hpp"
#include <math.h>
#include <vector>
#include <iostream>

class SphUtils {
public:
    //Default constructor
	SphUtils() {
    }

    //Explicit constructor
	SphUtils(float nsmooth_length, float nrho0, float nc, float nepsilon) {
		smooth_length = nsmooth_length;
		rho0 = nrho0;
		c = nc;
		epsilon = nepsilon;
    }

    //Destructor
	~SphUtils(){}

	//Kernel function for the SPH
	float kernel_function(glm::vec3 i, glm::vec3 j);
	glm::vec3 grad_kernel(glm::vec3 i, glm::vec3 j);

	void update_density(std::vector<Particle> &particles);
	void update_forces(std::vector<Particle> &particles);
	void update_posvel(std::vector<Particle> &particles, float dt);
	void collision(std::vector<Particle> &particles, std::vector<Wall> &walls);

private:
	float smooth_length;
	float rho0;
	float c;
	float epsilon;
	std::vector<std::pair<int,int>> collisions;

	void gravity_forces(std::vector<Particle> &particles);
	void pressure_forces(std::vector<Particle> &particles);
	void viscosity_forces(std::vector<Particle> &particles);
	void detect_collisions(std::vector<Particle> &particles, std::vector<Wall> &walls);
	void handle_collisions(std::vector<Particle> &particles, std::vector<Wall> &walls);
};