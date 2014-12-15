#define _USE_MATH_DEFINES

#include "SphUtils.h"
#include <glm/gtc/constants.hpp> 

/*
 * The kernel W(r,h). In this case, I used a cubic spline in 3 dimensions.
 */
float SphUtils::kernel_function(glm::vec3 i, glm::vec3 j) {
	//Distance in terms of smooth lenghts
	float q = glm::distance(i, j)/smooth_length;

	float result = -1.0f;
	if (q >= 2)      result = 0;
	else if (q > 1)  result = 1.0f/6.0f*(2-q)*(2-q)*(2-q);
	else if (q <= 1) result = 2.0f/3.0f - q*q + q*q*q/2;

	result *= 2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length);
	return result;
}

/*
 * The gradient of the kernel deltaW(r,h). In this case, I used a cubic spline in 3 dimensions.
 */
glm::vec3 SphUtils::grad_kernel(glm::vec3 i, glm::vec3 j) {
	//Returns the distance in terms of smooth lenghts
	float q = glm::distance(i, j)/smooth_length;
	float result = -1.0f;

	if (q >= 2)      result = 0;
	else if (q > 1)  result = -3.0f/6.0f*(2-q)*(2-q);
	else if (q <= 1) result = -2*q + 3*q*q/2;

	result *= 2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length);

	//Multiply the direction rij = ri - rj by the results
	glm::vec3 grad = glm::normalize(i-j);
	grad = grad* result;

	return grad;
}

/*
 * Updates the densities of all of the particles
 */
void SphUtils::update_density(std::vector<Particle> &particles) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		for (std::vector<Particle>::size_type j=0; j<particles.size(); j++) {
			//Pi = sum(mj * wij)
			particles[i].density += particles[j].mass*kernel_function(particles[i].pos, particles[j].pos);
		}
		//the pressure us updated using the speed of sound in water: pi = c^2 (rhoi - rho0)
		particles[i].pressure = c*c * (particles[i].density - rho0);
	}
}

/*
 * Updates the forces acting on all the particles.
 */
void SphUtils::update_forces(std::vector<Particle> &particles) {
	//Zero the forces
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		particles[i].force = glm::vec3(0.0f, 0.0f, 0.0f);
	}
	//Add the force from the pressure gradient
	pressure_forces(particles);
	//Add the gravity forces
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
 * The forces on the particles contributed by the pressures
 */
void SphUtils::pressure_forces(std::vector<Particle> &particles) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		for (std::vector<Particle>::size_type j=0; j<particles.size(); j++) {
			if ( i != j) {
				//deltaP/rho = F ~ -m*(pi/rhoi^2 + pj/rhoj^2)deltaW 
				float temp = -particles[j].mass *( particles[i].pressure/(particles[i].density*particles[i].density) + particles[j].pressure/(particles[j].density*particles[j].density) );
				particles[i].force += temp * grad_kernel(particles[i].pos, particles[j].pos);
			}
		}
	}
}

/*
 * Apply the gravitational force
 */
void SphUtils::gravity_forces(std::vector<Particle> &particles) {
	//std::cout << sph.kernel_function(particles[0].pos, particles[1].pos) << std::endl;
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//add the forces. For now, we only have simple gravity: F = mg
		particles[i].force += particles[i].mass*glm::vec3(0.0, -9.806, 0.0);
	}
}

/*
 * Detects and handles collisions
 */
void SphUtils::collision(std::vector<Particle> &particles, std::vector<Wall> &walls) {
	detect_collisions(particles,walls);
	handle_collisions(particles, walls);
}

/*
 * Detect collisions.
 */
void SphUtils::detect_collisions(std::vector<Particle> &particles, std::vector<Wall> &walls) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		for (std::vector<Wall>::size_type j=0; j<walls.size(); j++) {
			//Collision has occured if n dot (x - po) <= epsilon
			if (glm::dot(walls[j].normal, (particles[i].pos-walls[j].center)) <= epsilon) {
				std::cout << "Collision detected" << std::endl;
				collisions.push_back(std::make_pair(i,j));
			}
		}
	}
}

/*
 * Handle collisions. This is done by just applying an instantaneous impulse.
 */
void SphUtils::handle_collisions(std::vector<Particle> &particles, std::vector<Wall> &walls) {
	for (std::vector<std::pair<int,int>>::size_type i=0; i<collisions.size(); i++) {
		//Reflect the veloicty across the normal. Instan
		particles[collisions[i].first].vel = glm::reflect(particles[collisions[i].first].vel, walls[collisions[i].second].normal);
	}
}