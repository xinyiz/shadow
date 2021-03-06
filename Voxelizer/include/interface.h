#ifndef _GLUT_INTERFACE
#define _GLUT_INTERFACE
#define GL_GLEXT_PROTOTYPES
#include <GL/glut.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include "../include/vecmath/vecmath.h"
#include "../include/concQueue.h"

///////////////
// CONSTANTS //
///////////////
extern Vector3f addIndicator;
extern Vector3f remIndicator;
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
const int BRUSH_WIDTH = 6;

/////////////
// GLOBALS //
/////////////

// Lamp mesh data 
extern GLuint vbo_vertices;
extern GLuint vbo_normals;
extern GLuint vbo_colors;
extern GLuint ibo_elements;
extern GLuint ibo_elements_test;

// Lamp placement
extern float lamp_xpos;
extern float lamp_ypos;

// Cubic room data
extern GLuint room_vbo;
extern float room_dim;

// Light placement
extern float light_xpos;
extern float light_ypos;
extern float light_zpos;

typedef std::set<Vector3f> List;
// Drawn shadow pixels
// Sequence of point addition or removal actions. Each is a set of points to be 
// added or removed for the action. Addition actions are indicated by prescence 
// of (-1,-1) pair in set,  (-2, -2) for removal
extern Queue<List> shadowPixels;
// Tracking current pixel changes
extern List currShadowPixels;
// Points for drawing shadow in room 
extern std::set<Vector3f> pointsToDraw; 

// Interface parameters
extern int screenWidth;
extern int screenHeight;
extern bool mouseLeftDown;
extern bool mouseRightDown;
extern bool mouseDrawLeftDown;
extern bool mouseDrawRightDown;
extern int mouseX, mouseY;
extern int drawX, drawY;
extern float cameraAngleX;
extern float cameraAngleY;
extern float cameraDistance;
extern int brushWidth;
extern int mainWindow;
extern int drawWindow;
extern int sampleFactor;
extern int stride;

bool initSharedMem();
void initLights();
void setCamera(float posX, float posY, float posZ, float targetX, float targetY, float targetZ);
void toPerspective();
bool isValidPoint(int x, int y);
void displayCB();
void displayDrawCB();
int  initGLUT(int argc, char **argv);
int  initGLUTDraw(int argc, char **argv);
void initGL();
void timerCB(int millisec);
void timerDrawCB(int millisec);
void reshapeCB(int w, int h);
void mouseCB(int button, int stat, int x, int y);
Vector3f convertTo3DPoint(int x, int y);
void mouseDrawCB(int button, int stat, int x, int y);
void mouseMotionCB(int x, int y);
void mouseDrawMotionCB(int x, int y);
void keyboardCB(unsigned char key, int x, int y);
void deleteVBO(const GLuint vboId);
void clearSharedMem();
void exitCB();
#endif /* _GLUT_INTERFACE */
