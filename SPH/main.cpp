#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>

#include <GL/glew.h>
//#include "Particle.hpp"
#include <GL/glut.h>
//#include <GL/freeglut_ext.h>
//#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#define WINDOW_WIDTH 500
#define WINDOW_HEIGHT 500

GLuint programID;
GLuint vao;
GLuint vbo;

//std::vector<Particle> particles;

/*
 * Update function for GLUT
 *
 */	
void disp(void) {
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //particles[0].display();

	//glEnableVertexAttribArray(1);
	//Colors are indices [N...2N-1] in the vbo
	//glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*) ( WINDOW_WIDTH * WINDOW_HEIGHT * sizeof(glm::vec3)) );

	//glBindVertexArray(vao);
	//glDrawArrays(GL_POINTS, 0, 100);

	//unbind everything
	//glDisableVertexAttribArray(1);
	//glDisableVertexAttribArray(0);
	//glUseProgram(0);
	//glBindBuffer(GL_ARRAY_BUFFER, 0);
	//glFlush();
	//glutSwapBuffers();
}

/*
 * I don't think I need this
 *
 */ 
static void Reshape(int width, int height) {
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

/* Idle function */
static void idle( void ) {
	//int timems = glutGet(GLUT_ELAPSED_TIME);
	//if (timems % 100 == 0)
    //    glutPostRedisplay( );
}
/*
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
    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            Particle part(glm::vec3(i/10.0f,j/10.0f,0.0f), glm::vec3(0.0f,0.0f,0.0f), glm::vec3(0.0f,0.0f,0.0f), 1.0f, 1.0f, 1.0f, 1.0f);
            particles.push_back(part);
        }
    }
}
void init() {
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glEnable(GL_PROGRAM_POINT_SIZE);

	//Create the VAO	
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	//Load the shaders
	programID = LoadShaders( "sph.vertexshader", "sph.fragmentshader" );

	unsigned int size = 100;

    initScene();

    GLfloat initpos[3*size]; 
    int j=0;
    for (int i=0; i<particles.size(); i++) {
        initpos[j] = particles[i].pos.x;
        initpos[j+1] = particles[i].pos.y;
        initpos[j+2] = particles[i].pos.z;

        std::cout << initpos[j] << ", " << initpos[j+1] << ", " << initpos[j+2] << std::endl; 
        j += 3; 
    }

	//Create vertex buffer object
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	//Initialize the vbo to 0, as it will be computed by the GPU
	glBufferData(GL_ARRAY_BUFFER, 3*size, initpos, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	time_t start_time;
	time (&start_time);

	//register the vbo with opengl
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glUseProgram(programID);

	//Vertices are indices [0...N-1] in the vbo
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	glEnableVertexAttribArray(0);
}

*/
/*
 * Main function
 *
 */ 
int main (int argc, char** argv) {
	// init glut:
/*	glutInit (&argc, argv);
    glutInitContextVersion(3,3);
    glutInitContextProfile(GLUT_CORE_PROFILE);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition(0,0);
	glutCreateWindow("SPH");	
	//glutFullScreen();
    //glutGameModeString("1280x720:16@60"); glutEnterGameMode();
	//printf("OpenGL Version:%s\n",glGetString(GL_VERSION));
	//printf("GLSL Version  :%s\n",glGetString(GL_SHADING_LANGUAGE_VERSION));

	glutDisplayFunc(&disp);
    glutIdleFunc(&idle);
    glutKeyboardFunc(&keyboard);
	glutSpecialFunc(&skeyboard);
    glutMouseFunc(&mouse); 
*/
    // Initialise GLFW
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); //We don't want the old OpenGL
     
    // Open a window and create its OpenGL context
    GLFWwindow* window; // (In the accompanying source code, this variable is global)
    window = glfwCreateWindow( 500,500, "SPH", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window.\n" );
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window); // Initialize GLEW 

	glewExperimental=GL_TRUE;
	GLenum err=glewInit();
	if(err!=GLEW_OK)
		printf("glewInit failed, aborting.\n");
    if (!glewIsSupported("GL_VERSION_3_3 ")) {
        fprintf(stderr, "ERROR: Support for necessary OpenGL extensions missing.");
        fflush(stderr);
		exit(0);
    }

    std::cout << glGetString(GL_SHADING_LANGUAGE_VERSION ) << std::endl;
	//init();

    // Ensure we can capture the escape key being pressed below
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);     
    while (!glfwWindowShouldClose(window))
    {
        // Draw nothing, see you in tutorial 2 !
     
        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();
     
    } // Check if the ESC key was pressed or the window was closed

    glfwDestroyWindow(window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
	// enter tha main loop and process events:
//	glutMainLoop();
	return 0;
}
