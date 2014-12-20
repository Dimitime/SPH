#include "SphUtils.hpp"
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
	float result = 0.0f;

//	if (q >= 2)      result = 0;
//	else if (q > 1)  result = 1.0f/6.0f*(2-q)*(2-q)*(2-q);
//	else if (q <= 1) result = 2.0f/3.0f - q*q + q*q*q/2;
//	result *= 3.0f/(2.0f*(float)M_PI*smooth_length*smooth_length*smooth_length);

	//if (result != 0) std::cout << "Result nonzero: " << "q=" << q << " Wij= " << result*2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length) << std::endl;
	if (q <= smooth_length)
		result = 315* pow(smooth_length*smooth_length-q*q,3)/(64*(float)glm::pi<float>()*pow(smooth_length,9)) ;

	//std::cout << "r1: " << result << std::endl;
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

//	if (q >= 2)      result = 0;
///	else if (q > 1)  result = -3.0f/6.0f*(2-q)*(2-q);
//	else if (q <= 1) result = -2.0f*q + 3.0f*q*q/2.0f;
//	result *= 3.0f/(2.0f*(float)M_PI*smooth_length*smooth_length*smooth_length*smooth_length);

	if (q < smooth_length)
		result = -45.0f*(smooth_length-q)*(smooth_length-q)/((float)glm::pi<float>()*pow(smooth_length,6));

	//std::cout << result << std::endl;

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
/*	if (q >= 2)      result = 0;
	else if (q > 1)  result = -3.0f/6.0f*(2.0f-q)*(2.0f-q)-(2.0f-q)*r.x +
							  -3.0f/6.0f*(2.0f-q)*(2.0f-q)-(2.0f-q)*r.y + 
							  -3.0f/6.0f*(2.0f-q)*(2.0f-q)-(2.0f-q)*r.z;
	else if (q <= 1) result = -2.0f*q + 3.0f*q*q/2.0f+(-2.0f+3.0f*q)*r.x +
							  -2.0f*q + 3.0f*q*q/2.0f+(-2.0f+3.0f*q)*r.y +
							  -2.0f*q + 3.0f*q*q/2.0f+(-2.0f+3.0f*q)*r.z;
*/
	//if (result != 0) std::cout << "Result nonzero: " << "q=" << q << " Wij= " << result*2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length) << std::endl;
	if (q < smooth_length)
		result = 45/((float)glm::pi<float>()*pow(smooth_length,6)) * (smooth_length-q);
	
	//std::cout << "r1: " << result << std::endl;
	//result *= 2.0f/(3.0f*(float)M_PI*smooth_length*smooth_length*smooth_length*smooth_length);
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
//	for (std::vector<std::vector<int>>::size_type i=0; i<cells.size(); i++) {
		//loop through each particle in the cell
//		for (std::vector<int>::size_type m=0; m<cells[i].size(); m++) {
			
			
	//std::cout << "Min/Max: { (" << minx << ", " << maxx << ") " << " (" << miny << ", " << maxy << ")  " << " (" << minz << ", " << maxz << ")  " << std::endl;
	//std::cout << "Dimensions" << cell_dimensions << std::endl;
			for (size_t i=0; i<particles.size(); i++) {

			float density = 0.0f;

			int p1_index = i;//cells[i][m];
			for (size_t j=0; j<particles.size(); j++) {
			//We need to also check adjacent cells
/*			for (int a=-1; a<1; a++) {
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
*/								int p2_index = j;//cells[j][n];
								if (p1_index != p2_index) {			
									density += particles[p2_index].mass* glm::dot( (particles[p1_index].vel-particles[p2_index].vel), grad_kernel(particles[p1_index].pos, particles[p2_index].pos) );
													//std::cout << "Next!" << std::endl;
								}
//							}
//						}

					}
//				}
//			}
			//}
			particles[p1_index].density += dt*density;
			//particles[p1_index].pressure = 0.08f*(pow(particles[p1_index].density/rho0, 7)-1); //if (pow(particles[p1_index].density/rho0, 7)-1 > 2) std::cout << pow(particles[p1_index].density/rho0, 7)-1 << std::endl;
			//if (particles[p1_index].pressure < 0.0f)
			//	particles[p1_index].pressure = 0.0f;
			//the pressure us updated using the speed of sound in water: pi = c^2 (rhoi - rho0)
			particles[p1_index].pressure = 100.0f * (particles[p1_index].density - rho0);// std::cout << 10 * (particles[p1_index].density - rho0) << std::endl;
			//if (particles[p1_index].pressure > 100)
			//	std::cout << "Density for " << p1_index << ": " << particles[p1_index].density << " Pressure: " << particles[p1_index].pressure << std::endl;
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
//	for (std::vector<std::vector<int>>::size_type i=0; i<cells.size(); i++) {
//		//loop through each particle in the cell
//		for (std::vector<int>::size_type m=0; m<cells[i].size(); m++) {
			//and loop through the cells again...
			
			for (size_t i=0; i<particles.size(); i++) {
			int p1_index = i;//cells[i][m];
			for (size_t j=0; j<particles.size(); j++) {
			//We need to also check adjacent cells
/*			for (int a=-1; a<1; a++) {
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
*/								int p2_index = j;//cells[j][n];

								//std::cout << "Pairs" << i << ", " << j << " p1: " << p1_index << " p2: " << p2_index << std::endl;
								if ( p1_index != p2_index) {
									//std::cout << "BEFORE ACCESSING PARTICLE ARRAY" << std::endl;
									//deltaP/rho = F ~ -m*(pi/rhoi^2 + pj/rhoj^2)deltaW 
									//float temp = -particles[p2_index].mass *( particles[p1_index].pressure/(particles[p1_index].density*particles[p1_index].density) + particles[p2_index].pressure/(particles[p2_index].density*particles[p2_index].density) );

									//if (particles[p1_index].pressure*temp*glm::length(grad_kernel(particles[p1_index].pos, particles[p2_index].pos)) > 1)
									//	std::cout << /*particles[p1_index].pressure*/temp*glm::length(grad_kernel(particles[p1_index].pos, particles[p2_index].pos)) << std::endl;
									
									//We're going to add an artificial viscosity here
									float pi= 0.0f;// = // * grad_kernel(particles[p1_index].pos, particles[p2_index].pos)
									float vdotr = glm::dot( (particles[p1_index].vel-particles[p2_index].vel), (particles[p1_index].pos-particles[p2_index].pos));

									if (vdotr < 0) {
										float c = (sqrt(1.0f*abs(particles[p1_index].pressure/particles[p1_index].density))+sqrt(1.0f*abs(particles[p2_index].pressure/particles[p2_index].density)) )/2.0f;
										float r = glm::length(particles[p1_index].pos-particles[p2_index].pos);
										float mu = smooth_length* vdotr /( (r*r) + 0.01f*smooth_length*smooth_length);
										pi = 2.0f*(-0.1f*c*mu+0.2f*mu*mu)/ (particles[p1_index].density + particles[p2_index].density);
										//if (pi > 0.1)
										//std::cout << "Adding artifical viscosity of magintude " << pi << std::endl;
									}
									float temp = -particles[p2_index].mass*( pi + ( particles[p1_index].pressure+ particles[p2_index].pressure)/(2.0f*particles[p2_index].density*particles[p1_index].density));

									particles[p1_index].accel += temp * grad_kernel(particles[p1_index].pos, particles[p2_index].pos);
									//std::cout << "AFTER ACCESSING PARTICLE ARRAY" << std::endl;
								}
//							}
//						}
//					}
//				}
//			}
			//if (glm::length(particles[p1_index].accel) > 10)
			//	std::cout << "accel from pressure on " << p1_index << ": " << particles[p1_index].accel.x << ", " << particles[p1_index].accel.y << ", " << particles[p1_index].accel.z << std::endl;
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
//	for (std::vector<std::vector<int>>::size_type i=0; i<cells.size(); i++) {
		//loop through each particle in the cell
//		for (std::vector<int>::size_type m=0; m<cells[i].size(); m++) {
			//and loop through the cells again...
	for (size_t i=0; i<particles.size(); i++) {
			int p1_index = i;//cells[i][m];
	for (size_t j=0; j<particles.size(); j++) {

			//We need to also check adjacent cells
/*			for (int a=-1; a<1; a++) {
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
*/								int p2_index = j;//cells[j][n];

								//std::cout << "Pairs" << i << ", " << j << " p1: " << p1_index << " p2: " << p2_index << std::endl;
								if ( p1_index != p2_index) {
									//std::cout << "BEFORE ACCESSING PARTICLE ARRAY" << std::endl;
									//deltaP/rho = F ~ -m*(pi/rhoi^2 + pj/rhoj^2)deltaW 
									float temp = -0.1f*particles[p2_index].mass/(particles[p1_index].density*particles[p2_index].density) *lap_kernel(particles[p1_index].pos, particles[p2_index].pos);
									//std::cout << temp << "viscosity: " << std::endl;
									//if (temp*glm::length(particles[p1_index].vel-particles[p2_index].vel) > 1)
									//std::cout << "Force from visc: " << temp*glm::length(particles[p1_index].vel-particles[p2_index].vel) << std::endl;
									particles[p1_index].accel += temp*(particles[p1_index].vel-particles[p2_index].vel);
									//std::cout << "AFTER ACCESSING PARTICLE ARRAY" << std::endl;
								}
//							}
//						}
//					}
//				}
//			}
			//if (glm::length(particles[p1_index].force) > 10)
			//	std::cout << "Force from pressure on" << p1_index << ": " << particles[p1_index].force.x << ", " << particles[p1_index].force.y << ", " << particles[p1_index].force.z << std::endl;
		}
	}

}

/*
 * Apply the gravitational force
 */
void SphUtils::gravity_forces(std::vector<Particle> &particles,std::vector<Sphere> &spheres) {
	//std::cout << sph.kernel_function(particles[0].pos, particles[1].pos) << std::endl;
	for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
		//add the forces. For now, we only have simple gravity: F = mg
		particles[i].accel += /*particles[i].mass*/glm::vec3(0.0, -1, 0.0);

		//std::cout << "Total force on particle" << i << particles[i].force.x << ", " << particles[i].force.y << ", " <<particles[i].force.z << std::endl;
	}
	for (std::vector<Sphere>::size_type i=0; i<spheres.size(); i++) {
		//add the forces. For now, we only have simple gravity: F = mg
		spheres[i].accel += /*particles[i].mass*/glm::vec3(0.0, -1, 0.0);

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
		//Hack to keep it from exploding
		//if ( glm::dot(particles[i].vel,particles[i].vel) > 4) {
		//	particles[i].vel *= 4.0 / glm::dot(particles[i].vel,particles[i].vel);}
		//else {

		particles[i].vel = particles[i].vel +  dt * particles[i].accel;///particles[i].mass;
		//}
		//Use the updated velocities to calculate the new positions
		//xi+1 = xi + dt*xi+1

		//We want to adapt this for incompressable flow

		glm::vec3 vi = velocity_incompressible(i, particles);

		particles[i].pos = particles[i].pos +dt * (particles[i].vel);//+vi);
	}
	for (std::vector<Sphere>::size_type i=0; i<spheres.size(); i++) {
		//Update the velocities based on the forces
		// a = F/m
		//vi+1 = vi + dt*a
		//Hack to keep it from exploding
		//if ( glm::dot(particles[i].vel,particles[i].vel) > 4) {
		//	particles[i].vel *= 4.0 / glm::dot(particles[i].vel,particles[i].vel);}
		//else {

		spheres[i].vel = spheres[i].vel +  dt * spheres[i].accel;///particles[i].mass;
		//}
		//Use the updated velocities to calculate the new positions
		//xi+1 = xi + dt*xi+1

		//We want to adapt this for incompressable flow

		spheres[i].pos = spheres[i].pos +dt * (spheres[i].vel);//+vi);
	}
}

glm::vec3  SphUtils::velocity_incompressible(std::vector<Particle>::size_type i, std::vector<Particle> &particles) {
		glm::vec3 v = particles[i].vel;
		glm::vec3 vr(0.0f,0.0f,0.0f);
		float e = 0.5;
	
		for (size_t j=0; j<particles.size(); j++) {
			if (i != j)
			vr += 2*e*particles[j].mass * (particles[i].vel-particles[j].vel) * kernel_function(particles[i].pos,particles[j].pos) /(particles[i].density+particles[j].density);
			//std::cout << vr.length() << std::endl;
		}
		return vr;
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
			if ( dot+0.001 <= 0.1 || dot-0.001 <= 0.1 ) {
					//Reflect the veloicty across the normal.
					glm::vec3 impulse = glm::normalize(walls[j].normal);
					impulse *= (1.0f+0.0f)*glm::dot(walls[j].normal, particles[i].vel);

					//Repulsive force
					float fc = 10.0f*dot*glm::dot(-walls[j].normal, particles[i].accel);
					//std::cout << fc << std::endl;
					//std::cout << "Before " <<  particles[i].vel.x << ", " <<particles[i].vel.y << ", " <<particles[i].vel.z <<std::endl;

					if (fc < 0)
						fc = 0.0;
		
					glm::vec3 ac = fc*glm::normalize(walls[j].normal);
					//impulse *= 0.99;
					particles[i].vel += -impulse;
					//particles[i].accel += ac;
					
					//std::cout << "After " <<  particles[i].vel.x << ", " <<particles[i].vel.y << ", " <<particles[i].vel.z <<std::endl;
					//collisions.push_back(std::make_pair(i,j));
				//}
			}
		}

		//Check for collisions with there sphere
		for (std::vector<Sphere>::size_type j=0; j<spheres.size(); j++) {
			//Collision has occured if n dot (x - po) <= epsilon
			/*float dot = glm::dot(walls[j].normal, (particles[i].pos-walls[j].center));
			if ( dot+0.001 <= 0.1 || dot-0.001 <= 0.1 ) {
					//Reflect the veloicty across the normal.
					glm::vec3 impulse = glm::normalize(walls[j].normal);
					impulse *= (1.0f+0.0f)*glm::dot(walls[j].normal, particles[i].vel);

					//Repulsive force
					float fc = 10.0f*dot*glm::dot(-walls[j].normal, particles[i].accel);
					//std::cout << fc << std::endl;
					//std::cout << "Before " <<  particles[i].vel.x << ", " <<particles[i].vel.y << ", " <<particles[i].vel.z <<std::endl;

					if (fc < 0)
						fc = 0.0;
		
					glm::vec3 ac = fc*glm::normalize(walls[j].normal);
					//impulse *= 0.99;
					particles[i].vel += -impulse;
					//particles[i].accel += ac;
					
					//std::cout << "After " <<  particles[i].vel.x << ", " <<particles[i].vel.y << ", " <<particles[i].vel.z <<std::endl;
					//collisions.push_back(std::make_pair(i,j));
				//}
			}*/
			float dist = glm::length(particles[i].pos-spheres[j].pos);
			if (dist <= 0.1 + spheres[j].radius) {
				glm::vec3 u1 = particles[i].vel;
				glm::vec3 u2 = spheres[j].vel;
				particles[i].vel = (u1*(particles[i].mass-spheres[j].mass)+2.0f*spheres[j].mass*u2)/ (particles[i].mass+spheres[j].mass);
				std::cout << "New particle v: " << particles[i].vel.x <<  particles[i].vel.y << particles[i].vel.z << std::endl;
				spheres[j].vel = (u2*(spheres[j].mass-particles[i].mass)+2.0f*particles[i].mass*u1)/ (particles[i].mass+spheres[j].mass);
				std::cout << "New sphere v: " << spheres[j].vel.x <<  spheres[j].vel.y << spheres[j].vel.z << std::endl;
			}
		}
	}
}