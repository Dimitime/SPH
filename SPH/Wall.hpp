/*
 * Class definition for a Wall
 * the walls are pretty straightforward, but one thing to note is that in the initialization the pairs of
 * unconnected corners are c1/c4 and c2/c3
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

	void render_wall() {
	}

//private:
    glm::vec3 normal;
    glm::vec3 center;
	float xlength;
	float ylength;
};