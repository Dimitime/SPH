#include "SphUtils.hpp"
#include <glm/gtc/constants.hpp> 


/*
 * Build the nearest neighbors list
 */
void SphUtils::update_cells(std::vector<Particle> &particles) {
	//clear our previous cells
	cells.clear();
	maxx=minx=particles[0].pos.x;
	maxy=miny=particles[0].pos.y;
	maxz=minz=particles[0].pos.z;
	cell_dimensions = 2*smooth_length;
	//Loop through the particles to get the min and max values to determine how many cells we need.
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		if (particles[i].pos.x < minx)  minx = particles[i].pos.x;
		if (particles[i].pos.y < miny)  miny = particles[i].pos.y;
		if (particles[i].pos.z < minz)  minz = particles[i].pos.z;
		if (particles[i].pos.x > maxx)  maxx = particles[i].pos.x;
		if (particles[i].pos.y > maxy)  maxy = particles[i].pos.y;
		if (particles[i].pos.z > maxz)  maxz = particles[i].pos.z;
	}

	int nx = (int)ceil( (maxx-minx)/cell_dimensions)+1;
	int ny = (int)ceil( (maxy-miny)/cell_dimensions)+1;
	int nz = (int)ceil( (maxz-minz)/cell_dimensions)+1;
	size_t n =  nx*ny*nz;

	//Preallocate the size of the vector
	for (std::vector<std::vector<int>>::size_type i=0; i<n; i++) {
		std::vector<int> temp;
		cells.push_back(temp);
	}

	for (std::vector<Particle>::size_type m=0; m<particles.size(); m++) {
		//Turn the x,y,z coordinates of the particles into a cell index
		int x = (int)floor((particles[m].pos.x-minx)/cell_dimensions);
		int y = (int)floor((particles[m].pos.y-miny)/cell_dimensions);
		int z = (int)floor((particles[m].pos.z-minz)/cell_dimensions);

		//Get the flattened array index
		unsigned int i = x + nx*y + (nx*ny + ny)*z;
		
		cells[i].push_back(m);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// Kernel Functions
////////////////////////////////////////////////////////////////////////////////////////////
/*
 * The kernel W(r,h). In this case, I used a cubic spline in 3 dimensions.
 */
float SphUtils::kernel_function(glm::vec3 i, glm::vec3 j) {
	//Distance in terms of smooth lenghts
	float q = glm::distance(i, j)/smooth_length;
	float result = 0.0f;

	if (q <= smooth_length)
		result = 315* pow(smooth_length*smooth_length-q*q,3)/(64*(float)glm::pi<float>()*pow(smooth_length,9));
	return result;
}

/*
 * The gradient of the kernel deltaW(r,h). In this case, I used a cubic spline in 3 dimensions.
 */
glm::vec3 SphUtils::grad_kernel(glm::vec3 i, glm::vec3 j) {
	//Returns the distance in terms of smooth lenghts
	float q = glm::distance(i, j)/smooth_length;
	float result = 0.0f;

	if (q < smooth_length)
		result = -45.0f*(smooth_length-q)*(smooth_length-q)/((float)glm::pi<float>()*pow(smooth_length,6));
	glm::vec3 grad = glm::normalize(i-j);
	grad = grad* result;

	return grad;
}

/*
 * The laplacian of the kernel
 */
float SphUtils::lap_kernel(glm::vec3 i, glm::vec3 j) {
	//Distance in terms of smooth lenghts
	float q = glm::distance(i, j)/smooth_length;
	glm::vec3 r = glm::normalize(i-j);

	float result = 0.0f;
	if (q < smooth_length)
		result = 45/((float)glm::pi<float>()*pow(smooth_length,6)) * (smooth_length-q);
		//result = 45/((float)glm::pi<float>()*pow(smooth_length,5)) * (1-q/smooth_length);
	return result;

}

/*
 * Updates the densities of all of the particles
 */
void SphUtils::update_density(std::vector<Particle> &particles) {
	int nx = (int)ceil( (maxx-minx)/cell_dimensions)+1;
	int ny = (int)ceil( (maxy-miny)/cell_dimensions)+1;
	int nz = (int)ceil( (maxz-minz)/cell_dimensions)+1;

	//loop through the cells
//	for (std::vector<std::vector<int>>::size_type i=0; i<cells.size(); i++) {
	for (std::vector<std::vector<int>>::size_type i=0; i<particles.size(); i++) {
		//loop through each particle in the cell
//		for (std::vector<int>::size_type m=0; m<cells[i].size(); m++) {
			float density = 0.0f;

			int p1_index = i;//cells[i][m];
/*			//We need to also check adjacent cells
			for (int a=-1; a<1; a++) {
				for (int b=-1; b<1; b++) {
					for (int c=-1; c<1; c++) {
						int x = (int)floor((particles[p1_index].pos.x-minx)/cell_dimensions);
						int y = (int)floor((particles[p1_index].pos.y-miny)/cell_dimensions);
						int z = (int)floor((particles[p1_index].pos.z-minz)/cell_dimensions);
						//we need to clip the grids so they are valid indices
						int tx = x+a;
						int ty = y+b;
						int tz = z+c;

						//If this grid is out of bounds we skip it
						if ( (tx<0) || (tx>nx-1) || (ty<0) || (ty>ny-1) || (tz<0) || (tz>nz-1) ) {}
						else {
							unsigned int j = tx + nx*ty+ (nx*ny + ny)*tz;
							for (std::vector<int>::size_type n=0; n<cells[j].size(); n++) {
*/							for (std::vector<int>::size_type j=0; j<particles.size(); j++) {
								int p2_index = j;//cells[j][n];
								if (p1_index != p2_index) {			
									density += particles[p2_index].mass* kernel_function(particles[p1_index].pos, particles[p2_index].pos);//glm::dot( (particles[p1_index].vel-particles[p2_index].vel), grad_kernel(particles[p1_index].pos, particles[p2_index].pos) );
								}
							}
//						}

//					}
//				}
//			}
			if (density < rho0)
				density = rho0;
			particles[p1_index].density = density;
			particles[p1_index].pressure = 100.0f * (particles[p1_index].density - rho0);
			//if (p1_index == 1)
			//	std::cout << (particles[p1_index].density-rho0) << ", " << particles[p1_index].pressure << std::endl;
//		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// Force Functions
////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Updates the forces acting on all the particles.
 */
void SphUtils::update_forces(std::vector<Particle> &particles,std::vector<Sphere> &spheres) {
	//Zero the forces
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//add the forces. For now, we only have simple gravity: F = mg
			particles[i].accel = glm::vec3(0.0f, 0.0f, 0.0f);
	}
	for (std::vector<Sphere>::size_type i=0; i<spheres.size(); i++) {
		//add the forces. For now, we only have simple gravity: F = mg
		spheres[i].accel = glm::vec3(0.0f, 0.0f, 0.0f);
	}
	//Add the gravity forces
	gravity_forces(particles,spheres);
	//update_h_vel(particles);

	//Add the force from the pressure gradient
	pressure_forces(particles);

	//Add the forces due to viscosity
//	viscosity_forces(particles);
}

/*
 * The forces on the particles contributed by the pressures
 */
void SphUtils::pressure_forces(std::vector<Particle> &particles) {
	int nx = (int)ceil( (maxx-minx)/cell_dimensions)+1;
	int ny = (int)ceil( (maxy-miny)/cell_dimensions)+1;
	int nz = (int)ceil( (maxz-minz)/cell_dimensions)+1;

	//loop through the cells
//	for (std::vector<std::vector<int>>::size_type i=0; i<cells.size(); i++) {
	for (std::vector<std::vector<int>>::size_type i=0; i<particles.size(); i++) {
		//loop through each particle in the cell
//		for (std::vector<int>::size_type m=0; m<cells[i].size(); m++) {
			//and loop through the cells again...
			int p1_index = i;//cells[i][m];
/*			//We need to also check adjacent cells
			for (int a=-1; a<1; a++) {
				for (int b=-1; b<1; b++) {
					for (int c=-1; c<1; c++) {
						int x = (int)floor((particles[p1_index].pos.x-minx)/cell_dimensions);
						int y = (int)floor((particles[p1_index].pos.y-miny)/cell_dimensions);
						int z = (int)floor((particles[p1_index].pos.z-minz)/cell_dimensions);
						//we need to clip the grids so they are valid indices
						int tx = x+a;
						int ty = y+b;
						int tz = z+c;
						
						//If this grid is out of bounds we skip it
						if ( (tx<0) || (tx>nx-1) || (ty<0) || (ty>ny-1) || (tz<0) || (tz>nz-1) ) {}
						else {
							unsigned int j = tx + nx*ty+ (nx*ny + ny)*tz;

							//loop through each particle in the cell
							for (std::vector<int>::size_type n=0; n<cells[j].size(); n++) {
*/							for (std::vector<int>::size_type j=0; j<particles.size(); j++) {
								int p2_index = j;//cells[j][n];

								//std::cout << "Pairs" << i << ", " << j << " p1: " << p1_index << " p2: " << p2_index << std::endl;
								if ( p1_index != p2_index) {
									
									//We're going to add an artificial viscosity here
									float pi= 0.0f;// = // * grad_kernel(particles[p1_index].pos, particles[p2_index].pos)
									float vdotr = glm::dot( (particles[p1_index].vel-particles[p2_index].vel), (particles[p1_index].pos-particles[p2_index].pos));

									if (vdotr < 0) {
										float c = (sqrt(1.0f*abs(particles[p1_index].pressure/particles[p1_index].density))+sqrt(1.0f*abs(particles[p2_index].pressure/particles[p2_index].density)) )/2.0f;
										float r = glm::length(particles[p1_index].pos-particles[p2_index].pos);
										float mu = smooth_length* vdotr /( (r*r) + 0.01f*smooth_length*smooth_length);
										pi = 2.0f*(-0.1f*c*mu+0.2f*mu*mu)/ (particles[p1_index].density + particles[p2_index].density);
									}
									float temp = -particles[p2_index].mass*( pi + ( particles[p1_index].pressure+ particles[p2_index].pressure)/(2.0f*particles[p2_index].density*particles[p1_index].density));

									particles[p1_index].accel += temp * grad_kernel(particles[p1_index].pos, particles[p2_index].pos);
								}
							}
//						}
//					}
//				}
//			}
//		}
	}
}

/* 
 * The force from the viscosities of the particles
 */
void SphUtils::viscosity_forces(std::vector<Particle> &particles) {
	int nx = (int)ceil( (maxx-minx)/cell_dimensions)+1;
	int ny = (int)ceil( (maxy-miny)/cell_dimensions)+1;
	int nz = (int)ceil( (maxz-minz)/cell_dimensions)+1;

	//loop through the cells
	for (std::vector<std::vector<int>>::size_type i=0; i<cells.size(); i++) {
		//loop through each particle in the cell
		for (std::vector<int>::size_type m=0; m<cells[i].size(); m++) {
			//and loop through the cells again...
			int p1_index = cells[i][m];
			//We need to also check adjacent cells
			for (int a=-1; a<1; a++) {
				for (int b=-1; b<1; b++) {
					for (int c=-1; c<1; c++) {
						int x = (int)floor((particles[p1_index].pos.x-minx)/cell_dimensions);
						int y = (int)floor((particles[p1_index].pos.y-miny)/cell_dimensions);
						int z = (int)floor((particles[p1_index].pos.z-minz)/cell_dimensions);
						//we need to clip the grids so they are valid indices
						int tx = x+a;
						int ty = y+b;
						int tz = z+c;
						
						//If this grid is out of bounds we skip it
						if ( (tx<0) || (tx>nx-1) || (ty<0) || (ty>ny-1) || (tz<0) || (tz>nz-1) ) {}
						else {
							unsigned int j = tx + nx*ty+ (nx*ny + ny)*tz;

							//loop through each particle in the cell
							for (std::vector<int>::size_type n=0; n<cells[j].size(); n++) {
								int p2_index = cells[j][n];
								if ( p1_index != p2_index) {
									//std::cout << "BEFORE ACCESSING PARTICLE ARRAY" << std::endl;
									float temp = 0.1f*particles[p2_index].mass/(particles[p1_index].density*particles[p2_index].density) *lap_kernel(particles[p1_index].pos, particles[p2_index].pos);
									particles[p1_index].accel += temp*(particles[p1_index].vel-particles[p2_index].vel);
								}
							}
						}
					}
				}
			}
			//if (glm::length(particles[p1_index].force) > 10)
			//	std::cout << "Force from pressure on" << p1_index << ": " << particles[p1_index].force.x << ", " << particles[p1_index].force.y << ", " << particles[p1_index].force.z << std::endl;
		}
	}

}

/*
 * Apply the gravitational force
 */
void SphUtils::gravity_forces(std::vector<Particle> &particles,std::vector<Sphere> &spheres) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//add the forces. For now, we only have simple gravity
		particles[i].accel += /*particles[i].mass*/glm::vec3(0.0, -1, 0.0);
	}
	for (std::vector<Sphere>::size_type i=0; i<spheres.size(); i++) {
		//add the forces. For now, we only have simple gravity
		spheres[i].accel += /*particles[i].mass*/glm::vec3(0.0, -1, 0.0);
	}
}

/*
 * Updates the half timestep velocities for the leapfrog integration. Non-pressure forces are used.
 */
void SphUtils::update_h_vel(std::vector<Particle> &particles) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//Update the velocities based on the forces
		// a = F/m
		//vi+1 = vi + dt*a
		//particles[i].vel = particles[i].vel + dt/2 * particles[i].accel;///particles[i].mass;
	}
}

/*
 * Updates the positions and velocities of the particles
 */
void SphUtils::update_posvel(std::vector<Particle> &particles,std::vector<Sphere> &spheres) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//Update the velocities based on the forces
		// a = F/m
		//vi+1 = vi + dt*a
		particles[i].vel = particles[i].vel + dt*particles[i].accel;

		particles[i].pos = particles[i].pos + dt*(particles[i].vel);
		//std::cout << particles[i].pos.x << ","<< particles[i].pos.y << ","<< particles[i].pos.z << std::endl;
	}
	for (std::vector<Sphere>::size_type i=0; i<spheres.size(); i++) {
		spheres[i].vel = spheres[i].vel + dt*spheres[i].accel;
		spheres[i].pos = spheres[i].pos + dt*(spheres[i].vel);
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// Collision Functions
////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Detects and handles collisions
 */
void SphUtils::collision(std::vector<Particle> &particles, std::vector<Wall> &walls,std::vector<Sphere> &spheres) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//Check for collisions with the walls
		for (std::vector<Wall>::size_type j=0; j<walls.size(); j++) {
			//Collision has occured if n dot (x - po) <= epsilon
			float dot = glm::dot(walls[j].normal, (particles[i].pos-walls[j].center));
			if ( dot+0.001 <= 0.01 || dot-0.001 <= 0.01 ) {
					//Reflect the veloicty across the normal.
					glm::vec3 impulse = glm::normalize(walls[j].normal);
					impulse *= (1.0f+0.0f)*glm::dot(walls[j].normal, particles[i].vel);

					//Repulsive force
//					float fc = 10.0f*(0.01 - dot);//*glm::dot(-walls[j].normal, particles[i].accel);
//					if (fc < 0)
//						fc = 0.0;
//					
//					particles[i].accel += fc*walls[j].normal;

					//glm::vec3 ac = fc*glm::normalize(walls[j].normal);
					//impulse *= 0.99;

					//Put the particle back in legal range
					if (dot < 0)
						dot = -dot;
					particles[i].pos += dot*walls[j].normal;
					//std::cout << "pos: " << particles[i].pos.x << " " << particles[i].pos.y << " " << particles[i].pos.z << std::endl;
					particles[i].vel += -impulse;
					//std::cout << "vel: " << particles[i].vel.x << " " << particles[i].vel.y << " " << particles[i].vel.z << std::endl;
			}
		}

		//Check for collisions with there sphere
		for (std::vector<Sphere>::size_type j=0; j<spheres.size(); j++) {
			float dist = glm::length(particles[i].pos-spheres[j].pos);
			if (dist <= 0.1 + spheres[j].radius) {
				glm::vec3 u1 = particles[i].vel;
				glm::vec3 u2 = spheres[j].vel;
				particles[i].vel = (u1*(particles[i].mass-spheres[j].mass)+2.0f*spheres[j].mass*u2)/ (particles[i].mass+spheres[j].mass);
//				std::cout << "New particle v: " << particles[i].vel.x <<  particles[i].vel.y << particles[i].vel.z << std::endl;
				spheres[j].vel = (u2*(spheres[j].mass-particles[i].mass)+2.0f*particles[i].mass*u1)/ (particles[i].mass+spheres[j].mass);
//				std::cout << "New sphere v: " << spheres[j].vel.x <<  spheres[j].vel.y << spheres[j].vel.z << std::endl;
			}
		}
	}
}