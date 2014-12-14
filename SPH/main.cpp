#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <time.h>

#include <GL/glew.h>
#include <GL/glut.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>

#include "SphUtils.h"

#define WINDOW_WIDTH 500
#define WINDOW_HEIGHT 500
#define NUMBER_PARTICLES 100
#define NUMBER_WALLS 3

//OGL Buffer objects
GLuint programID;
GLuint vao;
GLuint vbo;
GLuint ebo;
GLuint MatrixID;

//timestep value
const float dt = 0.01f;
//Smoothing length
const float smooth_length = 0.1f;
//The ideal density. This is the density of water
const float rho0 = 1000.0f;
//The speed of sound in water
const float c = 10.0f;

//MVP matrices
// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
glm::mat4 Projection = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);

// Camera matrix
glm::mat4 View = glm::lookAt(
    glm::vec3(0,0,3), // Camera is at (4,3,3), in World Space
    glm::vec3(0,0,0), // and looks at the origin
    glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
);

// Model matrix : an identity matrix (model will be at the origin)
glm::mat4 Model = glm::mat4(1.0f);  // Changes for each model !

// Our ModelViewProjection : multiplication of our 3 matrices
glm::mat4 MVP = Projection * View * Model; // Remember, matrix multiplication is the other way around

SphUtils sph(smooth_length, rho0, c);
std::vector<Particle> particles;
std::vector<Wall> walls;

GLuint offset = 0;
GLuint elements[] = {
	offset,offset+1,
	offset+1,offset+2,
	offset+2,offset+3,
	offset+3,offset,

	offset+4,offset+5,
	offset+5,offset+6,
	offset+6,offset+7,
	offset+7,offset+4,

	offset+8,offset+9,
	offset+9,offset+10,
	offset+10,offset+11,
	offset+11,offset+8
};

/*
 * Draws the particles and the walls
 */
void draw_particles() {
	//Extract the positions so we can intialize the particles

	GLfloat* data;
	data = (GLfloat*) glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);
	if (data != (GLfloat*) NULL) {
		int j=0;
		for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
			data[j] = particles[i].pos.x;
			data[j+1] = particles[i].pos.y;
			data[j+2] = particles[i].pos.z;
			//std::cout << initpos[j] << ", " << initpos[j+1] << ", " << initpos[j+2] << std::endl; 
			j+=3; 
		}
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
}

void draw_walls() {

}

/*
 * Update function for GLUT
 *
 */	
void disp(void) {
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //particles[0].display();

	//glEnableVertexAttribArray(1);
	//Colors are indices [N...2N-1] in the vbo
	//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*) ( WINDOW_WIDTH * WINDOW_HEIGHT * sizeof(glm::vec3)) );
	
	glUseProgram(programID);
	
	//glBindBuffer(GL_ARRAY_BUFFER, vbo);	
	draw_particles();
	draw_walls();
 
	// Send our transformation to the currently bound shader,
	// in the "MVP" uniform
	// For each model you render, since the MVP will be different (at least the M part)
	glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

	//Vertices are indices [0...N-1] in the vbo
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(0);

	glBindVertexArray(vao);
	glDrawArrays(GL_POINTS, 0, NUMBER_PARTICLES);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glDrawElements(GL_POINTS,4*NUMBER_WALLS,GL_UNSIGNED_INT,elements);
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	//glDrawArrays(GL_LINES, NUMBER_PARTICLES, 4*NUMBER_WALLS);

	//unbind everything
	//glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(0);
	glUseProgram(0);
	//glBindBuffer(GL_ARRAY_BUFFER, 0);
	//glFlush();
	glutSwapBuffers();
}

////////////////////////////////////////////////////////////////////////////////////////////
// Control Functions
////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Mouse function for GLUT
 *
 */ 
void mouse(int button, int state, int x, int y) {
}

/*
 * Keyboard functions for GLUT
 *
 */ 
static void keyboard(unsigned char key, int x, int y) {
    switch (key) {
      case 27:
	exit(0);
    }
}

/*
 * Keyboard functions for GLUT
 *
 */ 
static void skeyboard(int key, int x, int y) {
}

///TODO: move to a different file
////////////////////////////////////////////////////////////////////////////////////////////
// Physics Functions
////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Advances the scene forward based on a symplectic euler integrator
 */
void step_scene() {
	std::cout << "Step" << std::endl;
	sph.update_density(particles);
	sph.update_forces(particles);
	sph.update_posvel(particles, dt);
}

/* Idle function */
static void idle() {
	int timems = glutGet(GLUT_ELAPSED_TIME);

	if (timems % 100 == 0) {
		step_scene();
		//std::cout << particles[0].pos.x << std::endl;
        glutPostRedisplay();
	}
}

//Loading shelders
GLuint LoadShaders(const char * vertex_file_path, const char * fragment_file_path) {
    // Create the shaders
    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);
 
    // Read the Vertex Shader code from the file
    std::string VertexShaderCode;
    std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
    if(VertexShaderStream.is_open())
    {
        std::string Line = "";
        while(std::getline(VertexShaderStream, Line))
            VertexShaderCode += "\n" + Line;
        VertexShaderStream.close();
    }
 
    // Read the Fragment Shader code from the file
    std::string FragmentShaderCode;
    std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
    if(FragmentShaderStream.is_open()){
        std::string Line = "";
        while(std::getline(FragmentShaderStream, Line))
            FragmentShaderCode += "\n" + Line;
        FragmentShaderStream.close();
    }
 
    GLint Result = GL_FALSE;
    int InfoLogLength;
 
    // Compile Vertex Shader
    printf("Compiling shader : %s\n", vertex_file_path);
    char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , NULL);
    glCompileShader(VertexShaderID);
 
    // Check Vertex Shader
    glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> VertexShaderErrorMessage(InfoLogLength);
    glGetShaderInfoLog(VertexShaderID, InfoLogLength, NULL, &VertexShaderErrorMessage[0]);
    fprintf(stdout, "%s\n", &VertexShaderErrorMessage[0]);
 
    // Compile Fragment Shader
    printf("Compiling shader : %s\n", fragment_file_path);
    char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , NULL);
    glCompileShader(FragmentShaderID);
 
    // Check Fragment Shader
    glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
    glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> FragmentShaderErrorMessage(InfoLogLength);
    glGetShaderInfoLog(FragmentShaderID, InfoLogLength, NULL, &FragmentShaderErrorMessage[0]);
    fprintf(stdout, "%s\n", &FragmentShaderErrorMessage[0]);
 
    // Link the program
    fprintf(stdout, "Linking program\n");
    GLuint ProgramID = glCreateProgram();
    glAttachShader(ProgramID, VertexShaderID);
    glAttachShader(ProgramID, FragmentShaderID);
    glLinkProgram(ProgramID);
 
    // Check the program
    glGetProgramiv(ProgramID, GL_LINK_STATUS, &Result);
    glGetProgramiv(ProgramID, GL_INFO_LOG_LENGTH, &InfoLogLength);
    std::vector<char> ProgramErrorMessage( glm::max(InfoLogLength, int(1)) );
    glGetProgramInfoLog(ProgramID, InfoLogLength, NULL, &ProgramErrorMessage[0]);
    fprintf(stdout, "%s\n", &ProgramErrorMessage[0]);
 
    glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);

    return ProgramID;
}

//Init the values ofthe particles
void initParticles() {
	//Density of water kg/m^3
	float density = rho0;
	//Start with 1 m^3 of water
	float volume = 1.0f;
	//Mass in KG of each particle
	float mass = density * volume / NUMBER_PARTICLES;
	
	//Pressure of the fluid
	float pressure = 1.0f;

	//TODO
	//Thermal energy? I might get rid of this
	float thermal = 1.0f;

    for (int i=-5; i<5; i++) {
        for (int j=-5; j<5; j++) {
			//for (int k=-5; k<5; k++) {
			float vel = -0.5f;
			if (i < 0)
				vel = 0.5f;
            Particle part(glm::vec3(i/10.0f,j/10.0f,0.0f), glm::vec3(0.0f,0.0f,0.0f), glm::vec3(0.0f,0.0f,0.0f), mass, density, pressure, thermal);
            particles.push_back(part);
			//}
        }
    }
}

//Init the values of the walls
void initWalls() {
	glm::vec3 c1(-0.6f,0.6f,0.6f);
	glm::vec3 c2(-0.6f,-1.3f,0.6f);
	glm::vec3 c3(-0.6f,0.6f,-0.6f);
	glm::vec3 c4(-0.6f,-1.3f,-0.6f);
	glm::vec3 c5(0.6f,0.6f,-0.60f);
	glm::vec3 c6(0.6f,-1.3f,-0.6f);
	glm::vec3 c7(0.6f,0.6f,0.6f);
	glm::vec3 c8(0.6f,-1.3f,0.6f);

	Wall w1(c1,c2,c3,c4);
	Wall w2(c4,c2,c6,c8);
	Wall w3(c5,c6,c7,c8);

	walls.push_back(w1);
	walls.push_back(w2);
	walls.push_back(w3);
}

//Initialize all the particles and walls
void initScene() {
	initParticles();
	initWalls();
}

//Initialize OpenGL
void init() {
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	//Create the VAO	
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	//Load the shaders
	programID = LoadShaders( "Shaders/sph.vertexshader", "Shaders/sph.fragmentshader" );

	// Get a handle for our "MVP" uniform.
	MatrixID = glGetUniformLocation(programID, "MVP");

    initScene();

	//Extract the positions so we can intialize the particles
    GLfloat initpos[3*(NUMBER_PARTICLES+4*NUMBER_WALLS)];
    int j=0;
    for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
        initpos[j] = particles[i].pos.x;
        initpos[j+1] = particles[i].pos.y;
        initpos[j+2] = particles[i].pos.z;
        //std::cout << initpos[j] << ", " << initpos[j+1] << ", " << initpos[j+2] << std::endl; 
        j += 3; 
    }

	std::cout << j << std::endl;
	for (std::vector<Wall>::size_type i=0; i<walls.size(); i++) {
		initpos[j] = walls[i].c1.x;
        initpos[j+1] = walls[i].c1.y;
        initpos[j+2] = walls[i].c1.z;
	
		initpos[j+3] = walls[i].c2.x;
        initpos[j+4] = walls[i].c2.y;
        initpos[j+5] = walls[i].c2.z;

		initpos[j+6] = walls[i].c3.x;
        initpos[j+7] = walls[i].c3.y;
        initpos[j+8] = walls[i].c3.z;

		initpos[j+9] = walls[i].c4.x;
        initpos[j+10] = walls[i].c4.y;
        initpos[j+11] = walls[i].c4.z;
/*
		initpos[j+12] = walls[i].c1.x;
        initpos[j+13] = walls[i].c1.y;
        initpos[j+14] = walls[i].c1.z;*/
        j += 12; 
	}
	std::cout << j << std::endl;

	//Create vertex buffer object
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	
	size_t size = 3*sizeof(GLfloat)*(NUMBER_PARTICLES+4*NUMBER_WALLS);
	//Initialize the vbo to the initial positions
	glBufferData(GL_ARRAY_BUFFER, size, initpos, GL_STREAM_DRAW);
	//glBindBuffer(GL_ARRAY_BUFFER, 0);

	int offset=NUMBER_PARTICLES;

	glGenBuffers(1, &ebo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);
}

/*
 * Main function
 *
 */ 
int main (int argc, char** argv) {
	// init glut:
	glutInit (&argc, argv);
    //glutInitContextVersion(3,3);
    //glutInitContextProfile(GLUT_CORE_PROFILE);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition(0,0);
	glutCreateWindow("SPH");	
	//glutFullScreen();
    //glutGameModeString("1280x720:16@60"); glutEnterGameMode();
	printf("OpenGL Version:%s\n",glGetString(GL_VERSION));
	printf("GLSL Version  :%s\n",glGetString(GL_SHADING_LANGUAGE_VERSION));

	glutDisplayFunc(&disp);
    glutIdleFunc(&idle);
    glutKeyboardFunc(&keyboard);
	glutSpecialFunc(&skeyboard);
    glutMouseFunc(&mouse); 

	glewExperimental=GL_TRUE;
	GLenum err=glewInit();
	if(err!=GLEW_OK)
		printf("glewInit failed, aborting.\n");
    if (!glewIsSupported("GL_VERSION_3_3 ")) {
        fprintf(stderr, "ERROR: Support for necessary OpenGL extensions missing.");
        fflush(stderr);
		exit(0);
    }

	init();

	// enter tha main loop and process events:
	glutMainLoop();
	return 0;
}