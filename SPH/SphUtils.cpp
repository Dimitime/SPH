#define _USE_MATH_DEFINES

#include "SphUtils.h"
#include <glm/gtc/constants.hpp> 


/*
 * Build the nearest neighbors list
 */
void SphUtils::update_cells(std::vector<Particle> &particles) {
	//clear our previous cells
	cells.clear();
	//float maxx, maxy, maxz, minx, miny, minz;
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

	//std::cout << "Min/Max: { (" << minx << ", " << maxx << ") " << " (" << miny << ", " << maxy << ")  " << " (" << minz << ", " << maxz << ")  " << std::endl;
	//std::cout << "Dimensions" << cell_dimensions << std::endl;

	int nx = (int)ceil( (maxx-minx)/cell_dimensions)+1;
	int ny = (int)ceil( (maxy-miny)/cell_dimensions)+1;
	int nz = (int)ceil( (maxz-minz)/cell_dimensions)+1;
	size_t n =  nx*ny*nz;

	//std::cout << "We need " << n << "=" << nx << "*" << ny << "*" << nz << " buckets" << std::endl;

	//Preallocate the size of the vector
	for (std::vector<std::vector<int>>::size_type i=0; i<n; i++) {
		std::vector<int> temp;
		cells.push_back(temp);
	}

	//std::cout << cells.size() << std::endl;

	for (std::vector<Particle>::size_type m=0; m<particles.size(); m++) {
		//Turn the x,y,z coordinates of the particles into a cell index
		int x = (int)floor((particles[m].pos.x-minx)/cell_dimensions);
		int y = (int)floor((particles[m].pos.y-miny)/cell_dimensions);
		int z = (int)floor((particles[m].pos.z-minz)/cell_dimensions);

		//Get the flattened array index
		unsigned int i = x + nx*y + (nx*ny + ny)*z;

		//std::cout << "[" << m << "]" << i  << " = " << x << ", " << y << ", " << z <<std::endl;

		if (particles[m].pos.z > 8.78) {
			//std::cout << "Pos: " << particles[m].pos.x << ", " << particles[m].pos.y << ", " << particles[m].pos.z << std::endl;
			//std::cout << "Extracted xyz: " << x << ", " << y << ", " << z << std::endl;

			//std::cout << m << ", " << i << "/" << cells.size() << std::endl;
		}
		
		cells[i].push_back(m);
	}
/*
	int test = 0;
	//A sanity check that we accounted for all of the particles
	for (int i=0; i< cells.size(); i++)
	{
		for (int j=0; j< cells[i].size(); j++)
		{
			test++;
		}
		if (cells[i].size() != 0)
		std::cout << "bucket size: " << cells[i].size() << std::endl;
	}
	
	std::cout << "TEST: " << test <<std::endl;
*/

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

	//std::cout << "Q: " << q << std::endl;

	float result = 0.0f;
	if (q >= 2)      result = 0;
	else if (q > 1)  result = 1.0f/6.0f*(2-q)*(2-q)*(2-q);
	else if (q <= 1) result = 2.0f/3.0f - q*q + q*q*q/2;

	//if (result != 0) std::cout << "Result nonzero: " << "q=" << q << " Wij= " << result*2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length) << std::endl;
	//if (q < smooth_length)
	//	result = 315/(64*(float)M_PI*pow(smooth_length,9)) * pow(smooth_length*smooth_length-q*q,3);

	//std::cout << "r1: " << result << std::endl;
	result *= 2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length);
	//std::cout << "r2: " << result << std::endl;
	return result;
}

/*
 * The gradient of the kernel deltaW(r,h). In this case, I used a cubic spline in 3 dimensions.
 */
glm::vec3 SphUtils::grad_kernel(glm::vec3 i, glm::vec3 j) {
	//Returns the distance in terms of smooth lenghts
	float q = glm::distance(i, j)/smooth_length;
	float result = 0.0f;

	if (q >= 2)      result = 0;
	else if (q > 1)  result = -3.0f/6.0f*(2-q)*(2-q);
	else if (q <= 1) result = -2.0f*q + 3.0f*q*q/2.0f;

	//if (q < smooth_length)
	//	result = -15.0f/((float)M_PI*pow(smooth_length,7))*(smooth_length-q)*(smooth_length-q);

	//std::cout << result << std::endl;
	result *= 2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length*smooth_length);

	//Multiply the direction rij = ri - rj by the results
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

	//std::cout << "Q: " << q << std::endl;

	float result = 0.0f;
	if (q >= 2)      result = 0;
	else if (q > 1)  result = -3.0f/6.0f*(2.0f-q)*(2.0f-q)-(2.0f-q)*r.x +
							  -3.0f/6.0f*(2.0f-q)*(2.0f-q)-(2.0f-q)*r.y + 
							  -3.0f/6.0f*(2.0f-q)*(2.0f-q)-(2.0f-q)*r.z;
	else if (q <= 1) result = -2.0f*q + 3.0f*q*q/2.0f+(-2.0f+3.0f*q)*r.x +
							  -2.0f*q + 3.0f*q*q/2.0f+(-2.0f+3.0f*q)*r.y +
							  -2.0f*q + 3.0f*q*q/2.0f+(-2.0f+3.0f*q)*r.z;

	//if (result != 0) std::cout << "Result nonzero: " << "q=" << q << " Wij= " << result*2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length) << std::endl;
	//if (q < smooth_length)
	//	result = 315/(64*(float)M_PI*pow(smooth_length,9)) * pow(smooth_length*smooth_length-q*q,3);

	//std::cout << "r1: " << result << std::endl;
	result *= 2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length);
	//std::cout << "r2: " << result << std::endl;
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
	for (std::vector<std::vector<int>>::size_type i=0; i<cells.size(); i++) {
		//loop through each particle in the cell
		for (std::vector<int>::size_type m=0; m<cells[i].size(); m++) {
			

	//std::cout << "Min/Max: { (" << minx << ", " << maxx << ") " << " (" << miny << ", " << maxy << ")  " << " (" << minz << ", " << maxz << ")  " << std::endl;
	//std::cout << "Dimensions" << cell_dimensions << std::endl;
			float density = 0.0f;

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
							//if (p1_index == 0) {
							//	std::cout << "number of particles in this bucket: " << cells[j].size() << std::endl;
							//	std::cout << "Pos: " << particles[p1_index].pos.x << ", " << particles[p1_index].pos.y << ", " << particles[p1_index].pos.z << std::endl;
							//}
							for (std::vector<int>::size_type n=0; n<cells[j].size(); n++) {
								int p2_index = cells[j][n];
								density += particles[p2_index].mass*kernel_function(particles[p1_index].pos, particles[p2_index].pos);
								//std::cout << "Next!" << std::endl;
							}
						}

					}
				}
			}
			//}
			particles[p1_index].density = density;
			particles[p1_index].pressure = 1.0f*(pow(particles[p1_index].density/rho0, 7)-1); //if (pow(particles[p1_index].density/rho0, 7)-1 > 2) std::cout << pow(particles[p1_index].density/rho0, 7)-1 << std::endl;
			//the pressure us updated using the speed of sound in water: pi = c^2 (rhoi - rho0)
			//particles[p1_index].pressure = 10 * (particles[p1_index].density - rho0); std::cout << 10 * (particles[p1_index].density - rho0) << std::endl;
			//if (particles[p1_index].pressure > 100)
			//	std::cout << "Density for " << p1_index << ": " << particles[p1_index].density << "Pressure: " << particles[p1_index].pressure << std::endl;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// Force Functions
////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Updates the forces acting on all the particles.
 */
void SphUtils::update_forces(std::vector<Particle> &particles) {
	//Zero the forces
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//add the forces. For now, we only have simple gravity: F = mg
			particles[i].accel = glm::vec3(0.0f, 0.0f, 0.0f);
	}
	//Add the gravity forces
	gravity_forces(particles);
	update_h_vel(particles);

	//Add the force from the pressure gradient
	pressure_forces(particles);

	//Add the forces due to viscosity
	//viscosity_forces(particles);
}

/*
 * The forces on the particles contributed by the pressures
 */
void SphUtils::pressure_forces(std::vector<Particle> &particles) {
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

								//std::cout << "Pairs" << i << ", " << j << " p1: " << p1_index << " p2: " << p2_index << std::endl;
								if ( p1_index != p2_index) {
									//std::cout << "BEFORE ACCESSING PARTICLE ARRAY" << std::endl;
									//deltaP/rho = F ~ -m*(pi/rhoi^2 + pj/rhoj^2)deltaW 
									//float temp = -particles[p2_index].mass *( particles[p1_index].pressure/(particles[p1_index].density*particles[p1_index].density) + particles[p2_index].pressure/(particles[p2_index].density*particles[p2_index].density) );
									float temp = -particles[p2_index].mass *( particles[p1_index].pressure+ particles[p2_index].pressure)/(2.0f*particles[p2_index].density);
									//if (particles[p1_index].pressure*temp*glm::length(grad_kernel(particles[p1_index].pos, particles[p2_index].pos)) > 1)
									//	std::cout << /*particles[p1_index].pressure*/temp*glm::length(grad_kernel(particles[p1_index].pos, particles[p2_index].pos)) << std::endl;

									//We're going to add an artificial viscosity here

									glm::vec3 av;// = // * grad_kernel(particles[p1_index].pos, particles[p2_index].pos)


									particles[p1_index].accel += 1/particles[p1_index].density*temp * grad_kernel(particles[p1_index].pos, particles[p2_index].pos);
									//std::cout << "AFTER ACCESSING PARTICLE ARRAY" << std::endl;
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

								//std::cout << "Pairs" << i << ", " << j << " p1: " << p1_index << " p2: " << p2_index << std::endl;
								if ( p1_index != p2_index) {
									//std::cout << "BEFORE ACCESSING PARTICLE ARRAY" << std::endl;
									//deltaP/rho = F ~ -m*(pi/rhoi^2 + pj/rhoj^2)deltaW 
									float temp = -.001f*particles[p2_index].mass/(particles[p2_index].density*particles[p2_index].density) *lap_kernel(particles[p1_index].pos, particles[p2_index].pos);
									//if (temp*glm::length(particles[p1_index].vel-particles[p2_index].vel) > 1)
									//std::cout << "Force from visc: " << temp*glm::length(particles[p1_index].vel-particles[p2_index].vel) << std::endl;
									particles[p1_index].accel += temp*(particles[p1_index].vel-particles[p2_index].vel);
									//std::cout << "AFTER ACCESSING PARTICLE ARRAY" << std::endl;
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
void SphUtils::gravity_forces(std::vector<Particle> &particles) {
	//std::cout << sph.kernel_function(particles[0].pos, particles[1].pos) << std::endl;
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//add the forces. For now, we only have simple gravity: F = mg
		particles[i].accel += /*particles[i].mass*/glm::vec3(0.0, -9.806f, 0.0);

		//std::cout << "Total force on particle" << i << particles[i].force.x << ", " << particles[i].force.y << ", " <<particles[i].force.z << std::endl;
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
		particles[i].vel = particles[i].vel + dt/2.0f * particles[i].accel;///particles[i].mass;
	}
}

/*
 * Updates the positions and velocities of the particles
 */
void SphUtils::update_posvel(std::vector<Particle> &particles) {
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//Update the velocities based on the forces
		// a = F/m
		//vi+1 = vi + dt*a
		particles[i].vel = particles[i].vel + dt/2 * particles[i].accel;///particles[i].mass;
		//Use the updated velocities to calculate the new positions
		//xi+1 = xi + dt*xi+1
		particles[i].pos = particles[i].pos + dt * particles[i].vel;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// Collision Functions
////////////////////////////////////////////////////////////////////////////////////////////
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
			float dot = glm::dot(walls[j].normal, (particles[i].pos-walls[j].center));

			if ( dot+0.001 <= epsilon || dot-0.001 <= epsilon ){
				
			    //float cos = dot / (walls[j].normal.length()* glm::length( (particles[i].pos-walls[j].center) ));
				//float sin = sqrt(1-cos*cos);

				//float xl = glm::length( (particles[i].pos-walls[j].center)*sin);
				//float l = T.length();
				//We have to check if we're out of bounds of the wall
				//if ( (xl-smooth_length >= walls[j].xlength) || (xl+smooth_length <= walls[j].xlength) ||
				//	 (xl-smooth_length >= walls[j].ylength) || (xl+smooth_length <= walls[j].ylength) ) {
				//if ((j = 4) || (j==5)) {
					//std::cout << i << "colldiing with back/front wall "<< j << std::endl;
					//std::cout << (particles[i].pos-walls[j].center).x << ", " << (particles[i].pos-walls[j].center).y << ", " << (particles[i].pos-walls[j].center).z << std::endl;
					///std::cout << "Particle position: " << (particles[i].pos).x << ", " << (particles[i].pos).y << ", " << (particles[i].pos).z << std::endl;
					//std::cout << "Wall center: " << (walls[j].center).x << ", " << (walls[j].center).y << ", " << (walls[j].center).z << std::endl;
					//std::cout << "Wall normal center: " << (walls[j].normal).x << ", " << (walls[j].normal).y << ", " << (walls[j].normal).z << std::endl;
					//std::cout << dot << std::endl << std::endl;
				//}

					collisions.push_back(std::make_pair(i,j));
				//}
			}
		}
	}
}

/*
 * Handle collisions. This is done by just applying an instantaneous impulse.
 */
void SphUtils::handle_collisions(std::vector<Particle> &particles, std::vector<Wall> &walls) {
	for (std::vector<std::pair<int,int>>::size_type i=0; i<collisions.size(); i++) {
		//Reflect the veloicty across the normal.
		glm::vec3 impulse = glm::normalize(walls[collisions[i].second].normal);
		impulse *= (1+0.9)*glm::dot(walls[collisions[i].second].normal, particles[collisions[i].first].vel);
		
		//impulse *= 0.99;
		particles[collisions[i].first].vel = particles[collisions[i].first].vel-impulse;;
	}
}