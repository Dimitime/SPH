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
        c1 = glm::vec3(0.0f,0.0f,0.0f);
        c2 = glm::vec3(0.0f,0.0f,0.0f);
        c3 = glm::vec3(0.0f,0.0f,0.0f);
        c4 = glm::vec3(0.0f,0.0f,0.0f);
    }

    //Explicit constructor
	Wall(glm::vec3 nc1, glm::vec3 nc2, glm::vec3 nc3, glm::vec3 nc4) {
        c1 = nc1;
		c2 = nc2;
		c3 = nc3;
		c4 = nc4;
    }

    //Destructor
	~Wall(){}

	void render_wall() {
	}

//private:
    glm::vec3 c1;
    glm::vec3 c2;
    glm::vec3 c3;
    glm::vec3 c4;
};