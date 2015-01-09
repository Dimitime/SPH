#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <time.h>

#include <GL/glew.h>
#include <GL/glut.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform2.hpp>
#include "SphUtils.hpp"

#define WINDOW_WIDTH 500
#define WINDOW_HEIGHT 500
#define NUMBER_PARTICLES 1300
#define NUMBER_WALLS 5
#define NUMBER_SPHERES 1
#define SPHERE_SLICES 40

//OGL Buffer objects
GLuint fluidShader;
GLuint sphereShader;
GLuint boxShader;
GLuint vao;
GLuint vbo[2];
//GLuint ebo;
GLuint MatrixID;

//timestep value
const float dt = 0.001f;
float total_time = 0.0f;

const float radius = 0.1f;
//Smoothing length
const float smooth_length = 6*radius;
//The ideal density. This is the density of water
const float rho0 = 10.0f;
//The speed of sound in water
const float c = 100.0f;
//An error value used in collision detection
const float epsilon = 0.1f;

bool toggleSim = true;

//MVP matrices
// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
glm::mat4 Projection = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);
// Camera matrix
glm::mat4 View = glm::lookAt(
    glm::vec3(1.5f, 1.0f,3.0f), // Camera is at (1.5,1,3), in World Space
    glm::vec3(0.0f, 0.0f, 0.0f), // and looks at the origin
    glm::vec3(0.0f,1.0f,0.0f)  // Head is up (set to 0,-1,0 to look upside-down)
);

// Model matrix : an identity matrix (model will be at the origin)
glm::mat4 Model = glm::mat4(1.0f);  // Changes for each model !

// Our ModelViewProjection : multiplication of our 3 matrices
glm::mat4 MVP = Projection * View * Model; // Remember, matrix multiplication is the other way around

SphUtils sph(smooth_length, rho0, c, dt, epsilon);
std::vector<Particle> particles;
std::vector<Wall> walls;
std::vector<Sphere> spheres;

/* 
 * Creates a triangle sphere mesh given a center anda radius
 */
std::vector<GLfloat> createSphereMesh(glm::vec3 center, float radius) {
	std::vector<GLfloat> v;

	for(int i = 0; i < SPHERE_SLICES; i++) {
		float theta0 = glm::pi<float>() * (-0.5f + (float) (i - 1) / SPHERE_SLICES);
		float theta1 = glm::pi<float>() * (-0.5f + (float) i / SPHERE_SLICES);

		for(int j = 0; j < SPHERE_SLICES; j++) {
			float phi = 2 * glm::pi<float>() * (float) (j - 1) / SPHERE_SLICES;
			v.push_back(radius*cos(phi)*cos(theta0)+center.x);
			v.push_back(radius*sin(phi)*cos(theta0)+center.y);
			v.push_back(radius* sin(theta0)+center.z);

			v.push_back(radius*cos(phi)*cos(theta1)+center.x);
			v.push_back(radius*sin(phi)*cos(theta1)+center.y);
			v.push_back(radius*sin(theta1)+center.z);
		}
	}

	return v;
}

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
	
		//Now get the sphere
		for (std::vector<Sphere>::size_type i=0; i<spheres.size(); i++) {
			//Construct a triangular mesh for every sphere
			std::vector<GLfloat> v = createSphereMesh(spheres[i].pos,spheres[i].radius);
			for (std::vector<GLfloat>::size_type k=0; k<v.size(); k+=3) {
				data[j] = v[k];//spheres[k].pos.x;
				data[j+1] = v[k+1]; //spheres[k].pos.y;
				data[j+2] = v[k+2];//spheres[k].pos.z;
				j+=3;
			}
		}
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
}

////////////////////////////////////////////////////////////////////////////////////////////
// Rendering and Animation Functions
////////////////////////////////////////////////////////////////////////////////////////////
/*
 * display function for GLUT
 */	
void disp(void) {
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	draw_particles();
	
	//Vertices are indices [0...N-1] in the vbo
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(0);
	//glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);	
	//Colors are indices [N...2N-1] in the vbo
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0 );
	glEnableVertexAttribArray(1);

	//Draw the fluid using the fluid shader
	glUseProgram(fluidShader);
	MatrixID = glGetUniformLocation(fluidShader, "MVP");
	glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	glDrawArrays(GL_POINTS, 0, NUMBER_PARTICLES);
	
	//Draw the spheres using the sphere shader
	//glUseProgram(sphereShader);
	//MatrixID = glGetUniformLocation(sphereShader, "MVP");
	//glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	glUseProgram(boxShader);
	MatrixID = glGetUniformLocation(boxShader, "MVP");
	glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	glDrawArrays(GL_TRIANGLE_STRIP, NUMBER_PARTICLES,2*SPHERE_SLICES*SPHERE_SLICES*NUMBER_SPHERES);
	
	//Draw the box using the box shader
	for (int i=0; i<NUMBER_WALLS; i++) {
		//std::cout << "drawing wall " << i << std::endl;
		glDrawArrays(GL_LINE_LOOP, NUMBER_PARTICLES+2*SPHERE_SLICES*SPHERE_SLICES*NUMBER_SPHERES+4*i, 4);
	}

	//unbind everything
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(0);
	glUseProgram(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glutSwapBuffers();
}

/*
 * Advances the scene forward based on a symplectic euler integrator
 */
void step_scene() {
	sph.update_cells(particles);
	//std::cout << "Grid updated" << std::endl;
	sph.update_density(particles);
	//std::cout << "Density updated" << std::endl;
	sph.update_forces(particles,spheres);
	//std::cout << "Forces updated" << std::endl;
	sph.update_posvel(particles,spheres);
	//std::cout << "Pos/Vel updated" << std::endl;
	sph.collision(particles, walls, spheres);

	//std::cout << "Step done" << std::endl;
	total_time += dt;
}

/* 
 * Idle function 
 */
static void idle() {
	int timems = glutGet(GLUT_ELAPSED_TIME);

	if (timems % 100 == 0) {
		if (toggleSim)
			step_scene();
        glutPostRedisplay();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////
// Control Functions
////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Keyboard functions for GLUT
 *
 */ 
static void keyboard(unsigned char key, int x, int y) {
    switch (key) {
		case 's':
			step_scene();
			glutPostRedisplay();
			break;
		case 32:
			toggleSim = !toggleSim;
			break;
		case 27:
	exit(0);
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
	float density = rho0+epsilon;
	//Mass in KG of each particle
	float mass = .1f*smooth_length*smooth_length*smooth_length*rho0;
	std::cout << "Mass: " << mass << std::endl;
	//Pressure of the fluid
	float pressure = 1.0f;

	//TODO
	//Thermal energy? I might get rid of this
	float thermal = 1.0f;

	float x = -0.5f;
	float zvel = 0.0f;
	float yvel = 0.0f;//-1.0f;
    for (int i=0; i<13; i++) {
		float y = -0.9f;
        for (int j=0; j<10; j++) {
			float z =  -0.5f;
			for (int k=0; k<10; k++) {
				Particle part(glm::vec3(x,y,z), glm::vec3(0.0f,yvel,zvel*j), glm::vec3(0.0f,0.0f,0.0f), mass, density, pressure, thermal);
				particles.push_back(part);
				z += radius;
			}
			y += radius;
        }
		x += radius;
    }

	std::cout << particles.size() << std::endl;
}

//Init the values of the walls
void initWalls() {
	//The side walls
	glm::vec3 center1(-1.0f, 0.0f, 0.0f);
	glm::vec3 normal1(1.0f, 0.0f, 0.0f);
	float xlength = 1.0f;
	float ylength = 1.0f;

	glm::vec3 center2(1.0f, 0.0f, 0.0f);
	glm::vec3 normal2(-1.0f, 0.0f, 0.0f);
	Wall w1(center1,normal1,xlength,ylength);
	Wall w2(center2,normal2,xlength,ylength);

	//The bottom wall
	glm::vec3 center3(0.0f, -1.0f, 0.0f);
	glm::vec3 normal3(0.0f, 1.0f, 0.0f);
	Wall w3(center3,normal3,xlength,xlength);

	//The front wall
	glm::vec3 center4(0.0f, 0.0f,1.0);
	glm::vec3 normal4(0.0f, 0.0f,-1.0f);
	Wall w4(center4,normal4,xlength,xlength);

	//The back wall
	glm::vec3 center5(0.0f, 0.0f, -1.0f);
	glm::vec3 normal5(0.0f, 0.0f,1.0f);
	Wall w5(center5,normal5,xlength,xlength);

	walls.push_back(w1);
	walls.push_back(w2);
	walls.push_back(w4);
	walls.push_back(w5);
	walls.push_back(w3);
}

void initSpheres() {
	glm::vec3 pos(0.0f, 0.25f, 0.0f);
	glm::vec3 vel(0.0f,0.0f, 0.0f);
	float raidus = 0.2f;
	float mass = 1.2f;

	Sphere s(pos, vel, radius, mass);
	spheres.push_back(s);
}

//Initialize all the particles and walls
void initScene() {
	initParticles();
	initSpheres();
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
	fluidShader = LoadShaders( "Shaders/sph.vertexshader", "Shaders/sph.fragmentshader" );
	boxShader = LoadShaders( "Shaders/sph.vertexshader", "Shaders/sph.fragmentshader" );
	//sphereShader = LoadShaders( "Shaders/sph.vertexshader", "Shaders/sph.fragmentshader" );

    initScene();
	GLuint sphereVerts = 2*SPHERE_SLICES*SPHERE_SLICES; //3200 vertixes
	//Extract the positions so we can intialize the particles
    GLfloat* initpos = new GLfloat[3*(NUMBER_PARTICLES+4*NUMBER_WALLS+sphereVerts*NUMBER_SPHERES)];
    int j=0;
    for (std::vector<Particle>::size_type i=0; i<particles.size(); i++) {
        initpos[j] = particles[i].pos.x;
        initpos[j+1] = particles[i].pos.y;
        initpos[j+2] = particles[i].pos.z; 
        j += 3; 
    }
	
	//Now get the sphere
	for (std::vector<Sphere>::size_type i=0; i<spheres.size(); i++) {
		//Construct a triangular mesh for every sphere
		std::vector<GLfloat> v = createSphereMesh(spheres[i].pos,spheres[i].radius);
		std::cout << "Number of vertices in sphere: " << v.size();
		for (std::vector<GLfloat>::size_type k=0; k<v.size(); k+=3) {
			initpos[j] = v[k];
			initpos[j+1] = v[k+1];
			initpos[j+2] = v[k+2];
			j+=3;
		}
	}

	//The two side walls
	for (std::vector<Wall>::size_type i=0; i<2; i++) {
		initpos[j] = walls[i].center.x;
        initpos[j+1] = walls[i].center.y+walls[i].ylength;
        initpos[j+2] = walls[i].center.z+walls[i].xlength;

		initpos[j+3] = walls[i].center.x;
        initpos[j+4] = walls[i].center.y-walls[i].ylength;
        initpos[j+5] = walls[i].center.z+walls[i].xlength;

		initpos[j+6] = walls[i].center.x;
        initpos[j+7] = walls[i].center.y-walls[i].ylength;
        initpos[j+8] = walls[i].center.z-walls[i].xlength;

		initpos[j+9] = walls[i].center.x;
        initpos[j+10] = walls[i].center.y+walls[i].ylength;
        initpos[j+11] = walls[i].center.z-walls[i].xlength;
        j += 12; 
	}

	//The front and back wall
	for (std::vector<Wall>::size_type i=2; i<4; i++) {
		initpos[j] = walls[i].center.x+walls[i].xlength;
        initpos[j+1] = walls[i].center.y+walls[i].ylength;
        initpos[j+2] = walls[i].center.z;
	
		initpos[j+3] = walls[i].center.x+walls[i].xlength;
        initpos[j+4] = walls[i].center.y-walls[i].ylength;
        initpos[j+5] = walls[i].center.z;

		initpos[j+6] = walls[i].center.x-walls[i].xlength;
        initpos[j+7] = walls[i].center.y-walls[i].ylength;
        initpos[j+8] = walls[i].center.z;

		initpos[j+9] = walls[i].center.x-walls[i].xlength;
        initpos[j+10] = walls[i].center.y+walls[i].ylength;
        initpos[j+11] = walls[i].center.z;
        j += 12; 
	}

	//The bottom wall
	int i = walls.size()-1;

	initpos[j] = walls[i].center.x+walls[i].xlength;
    initpos[j+1] = walls[i].center.y;
    initpos[j+2] = walls[i].center.z+walls[i].ylength;
	
	initpos[j+3] = walls[i].center.x+walls[i].xlength;
    initpos[j+4] = walls[i].center.y;
    initpos[j+5] = walls[i].center.z-walls[i].ylength;

	initpos[j+6] = walls[i].center.x-walls[i].xlength;
    initpos[j+7] = walls[i].center.y;
    initpos[j+8] = walls[i].center.z-walls[i].ylength;

	initpos[j+9] = walls[i].center.x-walls[i].xlength;
    initpos[j+10] = walls[i].center.y;
    initpos[j+11] = walls[i].center.z+walls[i].ylength;


	//Add colors
	j=0;
    GLfloat* colors = new GLfloat[3*(NUMBER_PARTICLES+4*NUMBER_WALLS+sphereVerts*NUMBER_SPHERES)];

	//Set the particles to be blue
	for (std::vector<Particle>::size_type c=0; c<particles.size(); c++) {
        colors[j] = 0;
		colors[j+1] = 0;
        colors[j+2] = 1.0f;
		j+=3;
	}

	//Set the sphere to be green.
	for (std::vector<Sphere>::size_type c=0; c<spheres.size(); c++) {
		//Construct a triangular mesh for every sphere
		std::vector<GLfloat> v = createSphereMesh(spheres[i].pos,spheres[i].radius);
		std::cout << "Number of vertices in sphere2: " << v.size();
		for (std::vector<GLfloat>::size_type k=0; k<v.size(); k+=3) {
			colors[j] = 0.0f;
			colors[j+1] = 1.0f;
			colors[j+2] = 0.0f;
		j+=3;
		}
	}

	//Set walls to be red?
	for (std::vector<Wall>::size_type c=0; c<4*walls.size(); c++) {
        colors[j] = 1.0f;
        colors[j+1] = 0;
        colors[j+2] = 0;
		j+=3;
	}

	//Create vertex buffer object
	glGenBuffers(2, vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
	
	size_t size1 = 3*sizeof(GLfloat)*(NUMBER_PARTICLES+4*NUMBER_WALLS+sphereVerts*NUMBER_SPHERES);
	//Initialize the vbo to the initial positions
	glBufferData(GL_ARRAY_BUFFER, size1, initpos, GL_STREAM_DRAW);

	//glGenBuffers(1, &vboc);
    glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
	size_t size2 = 3*sizeof(GLubyte)*(NUMBER_PARTICLES+4*NUMBER_WALLS);
	glBufferData(GL_ARRAY_BUFFER, size1, colors, GL_STREAM_DRAW);
    glEnableVertexAttribArray(1);

	delete [] initpos;
	delete [] colors;

	std::cout <<"Buffer setup complete" << std::endl;
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