#ifndef _GLUT_INTERFACE
#define _GLUT_INTERFACE
#define GL_GLEXT_PROTOTYPES
#include <GL/glut.h>

///////////////
// CONSTANTS //
///////////////
const int   SCREEN_WIDTH    = 400;
const int   SCREEN_HEIGHT   = 300;
const float CAMERA_DISTANCE = 5.0f;

/////////////
// GLOBALS //
/////////////

// Lamp Mesh Data 
extern GLuint vbo_vertices;
extern GLuint vbo_normals;
extern GLuint ibo_elements;

// Interface parameters
extern int screenWidth;
extern int screenHeight;
extern bool mouseLeftDown;
extern bool mouseRightDown;
extern float mouseX, mouseY;
extern float cameraAngleX;
extern float cameraAngleY;
extern float cameraDistance;

bool initSharedMem();
void initLights();
void setCamera(float posX, float posY, float posZ, float targetX, float targetY, float targetZ);
void toPerspective();
void displayCB();
int  initGLUT(int argc, char **argv);
void initGL();
void timerCB(int millisec);
void reshapeCB(int w, int h);
void mouseCB(int button, int stat, int x, int y);
void mouseMotionCB(int x, int y);
void deleteVBO(const GLuint vboId);
void clearSharedMem();
void exitCB();
#endif /* _GLUT_INTERFACE */
