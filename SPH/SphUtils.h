/*
 * Class definition for a Particle
 *
 */

#include "Particle.hpp"
#include <math.h>
#include <vector>
#include <iostream>

class SphUtils {
public:
    //Default constructor
	SphUtils() {
    }
/*
    //Explicit constructor
	SphUtils() {
    }
*/
    //Destructor
	~SphUtils(){}

	//Kernel function for the SPH
	float kernel_function(glm::vec3 i, glm::vec3 j, float smooth_length);

	void update_density(std::vector<Particle> &particles);
	void update_forces(std::vector<Particle> &particles);
	void update_posvel(std::vector<Particle> &particles, float dt);

private:
	void gravity_forces(std::vector<Particle> &particles);
};