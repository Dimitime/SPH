/*
 * Class definition for a Wall
 */

#include <glm/glm.hpp>

class Wall {
public:
    //Default constructor
	Wall() {
        center = glm::vec3(0.0f,0.0f,0.0f);
        normal = glm::vec3(0.0f,0.0f,0.0f);
        xlength = 0.0f;
		ylength = 0.0f;
    }

    //Explicit constructor
	Wall(glm::vec3 ncenter, glm::vec3 nnormal, float nxlength, float nylength) {
		center = ncenter;
		normal = nnormal;
		xlength = nxlength;
		ylength = nylength;
    }

    //Destructor
	~Wall(){}

//private:
    glm::vec3 normal;
    glm::vec3 center;
	float xlength;
	float ylength;
};