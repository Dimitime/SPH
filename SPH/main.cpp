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

//OGL Buffer objects
GLuint programID;
GLuint vao;
GLuint vbo;
GLuint MatrixID;

//timestep value
const float dt = 0.01f;
const float smooth_length = 1.0f;

//MVP matrices
// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
glm::mat4 Projection = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);

// Camera matrix
glm::mat4 View = glm::lookAt(
    glm::vec3(0,0,5), // Camera is at (4,3,3), in World Space
    glm::vec3(0,0,0), // and looks at the origin
    glm::vec3(0,1,0)  // Head is up (set to 0,-1,0 to look upside-down)
);

// Model matrix : an identity matrix (model will be at the origin)
glm::mat4 Model = glm::mat4(1.0f);  // Changes for each model !

// Our ModelViewProjection : multiplication of our 3 matrices
glm::mat4 MVP = Projection * View * Model; // Remember, matrix multiplication is the other way around

SphUtils sph;
std::vector<Particle> particles;

/*
 * Draws the particles
 */
void draw_particles() {
	//register the vbo with opengl
	glBindBuffer(GL_ARRAY_BUFFER, vbo);	
	//Extract the positions so we can intialize the particles
    GLfloat initpos[3*NUMBER_PARTICLES]; 
    int j=0;
    for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
        initpos[j] = particles[i].pos.x;
        initpos[j+1] = particles[i].pos.y;
        initpos[j+2] = particles[i].pos.z;
        //std::cout << initpos[j] << ", " << initpos[j+1] << ", " << initpos[j+2] << std::endl; 
        j += 3; 
    }
	//Initialize the vbo to 0, as it will be computed by the GPU
	glBufferData(GL_ARRAY_BUFFER, 3*NUMBER_PARTICLES*sizeof(glm::vec3), initpos, GL_DYNAMIC_DRAW);
}

/*
 * Update function for GLUT
 *
 */	
void disp(void) {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //particles[0].display();

	//glEnableVertexAttribArray(1);
	//Colors are indices [N...2N-1] in the vbo
	//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*) ( WINDOW_WIDTH * WINDOW_HEIGHT * sizeof(glm::vec3)) );
	
	draw_particles();

	glUseProgram(programID);
 
	// Send our transformation to the currently bound shader,
	// in the "MVP" uniform
	// For each model you render, since the MVP will be different (at least the M part)
	glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

	//Vertices are indices [0...N-1] in the vbo
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(0);

	glBindVertexArray(vao);
	glDrawArrays(GL_POINTS, 0, 100);

	//unbind everything
	//glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(0);
	glUseProgram(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
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

//Initialize all the particles
void initScene() {

	//Density of water kg/m^3
	float density = 999.97f;
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
            Particle part(glm::vec3(i/10.0f,j/10.0f,0.0f), glm::vec3(0.0f,0.0f,0.0f), glm::vec3(0.0f,0.0f,0.0f), mass, density, pressure, thermal);
            particles.push_back(part);
        }
    }
}

//Initialize OpenGL
void init() {
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glEnable(GL_PROGRAM_POINT_SIZE);

	//Create the VAO	
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	//Load the shaders
	programID = LoadShaders( "Shaders/sph.vertexshader", "Shaders/sph.fragmentshader" );

	// Get a handle for our "MVP" uniform.
	MatrixID = glGetUniformLocation(programID, "MVP");

	const unsigned int size = 100;

    initScene();

	//Extract the positions so we can intialize the particles
    GLfloat initpos[3*NUMBER_PARTICLES]; 
    int j=0;
    for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
        initpos[j] = particles[i].pos.x;
        initpos[j+1] = particles[i].pos.y;
        initpos[j+2] = particles[i].pos.z;
        //std::cout << initpos[j] << ", " << initpos[j+1] << ", " << initpos[j+2] << std::endl; 
        j += 3; 
    }
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	//Initialize the vbo to the initial positions
	glBufferData(GL_ARRAY_BUFFER, 3*NUMBER_PARTICLES*sizeof(glm::vec3), initpos, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//Create vertex buffer object
	glGenBuffers(1, &vbo);

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