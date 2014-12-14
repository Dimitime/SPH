#define _USE_MATH_DEFINES

#include "SphUtils.h"
#include <glm/gtc/constants.hpp> 

//A cubic spline kernel function
float SphUtils::kernel_function(glm::vec3 i, glm::vec3 j, float smooth_length) {
	//Returns the distance in terms of smooth lenghts
	float q = glm::distance(i, j)/smooth_length;
	float result = -1.0f;
	if (q >= 2)      result = 0;
	else if (q > 1)  result = 1.0f/6.0f*(2-q)*(2-q)*(2-q);
	else if (q <= 1) result = 2.0f/3.0f - q*q + q*q*q/2;

	result *= 2.0f/(3.0f*M_PI*smooth_length*smooth_length*smooth_length);
	return result;
}

/*
 * Updates the densities of all of the particles
 */
void SphUtils::update_density(std::vector<Particle> &particles) {
}

/*
 * Updates the forces acting on all the particles.
 */
void SphUtils::update_forces(std::vector<Particle> &particles) {
	gravity_forces(particles);
}

/*
 * Updates the positions of the particles
 */
void SphUtils::update_posvel(std::vector<Particle> &particles, float dt) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//Update the velocities based on the forces
		// a = F/m
		//vi+1 = vi + dt*a
		particles[i].vel = particles[i].vel + dt * particles[i].force/particles[i].mass;
		//Use the updated velocities to calculate the new positions
		//xi+1 = xi + dt*xi+1
		particles[i].pos = particles[i].pos + dt * particles[i].vel;
	}
}

/*
 * Apply the gravitational force
 */
void SphUtils::gravity_forces(std::vector<Particle> &particles) {
	//std::cout << sph.kernel_function(particles[0].pos, particles[1].pos) << std::endl;
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//add the forces. For now, we only have simple gravity: F = mg
		particles[i].force = particles[i].mass*glm::vec3(0.0, -9.806, 0.0);
	}
}