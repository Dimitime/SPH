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
	SphUtils(float nsmooth_length, float nrho0, float nc) {
		smooth_length = nsmooth_length;
		rho0 = nrho0;
		c = nc;
    }

    //Destructor
	~SphUtils(){}

	//Kernel function for the SPH
	float kernel_function(glm::vec3 i, glm::vec3 j);
	glm::vec3 grad_kernel(glm::vec3 i, glm::vec3 j);

	void update_density(std::vector<Particle> &particles);
	void update_forces(std::vector<Particle> &particles);
	void update_posvel(std::vector<Particle> &particles, float dt);

private:
	float smooth_length;
	float rho0;
	float c;
	void gravity_forces(std::vector<Particle> &particles);
	void pressure_forces(std::vector<Particle> &particles);
	void viscosity_forces(std::vector<Particle> &particles);
};