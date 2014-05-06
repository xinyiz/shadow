#include "../include/interface.h"

int screenWidth;
int screenHeight;
bool mouseLeftDown;
bool mouseRightDown;
bool mouseDrawLeftDown;
bool mouseDrawRightDown;
float mouseX, mouseY;
float drawX, drawY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance;
int brushWidth;
int mainWindow;
int drawWindow;
std::vector<List> shadowPixels;
List currShadowPixels;
vector3fList convertedShadowPixels;

///////////////////////////////////////////////////////////////////////////////
// initialize global variables
///////////////////////////////////////////////////////////////////////////////
bool initSharedMem()
{
    screenWidth = SCREEN_WIDTH;
    screenHeight = SCREEN_HEIGHT;
    mouseLeftDown = mouseRightDown = false;
    mouseX = mouseY = 0;
    drawX = drawY = 0;
    cameraAngleX = cameraAngleY = 0.0f;
    cameraDistance = room_dim*2.0f;
    brushWidth = BRUSH_WIDTH;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// initialize lights
///////////////////////////////////////////////////////////////////////////////
void initLights()
{
    // set up light colors (ambient, diffuse, specular)
    GLfloat lightKd[] = {0.2f, 0.6f, 0.9f, 1.0f};  // diffuse light
    GLfloat lightKs[] = {1, 1, 1, 1};           // specular light
    //glLightfv(GL_LIGHT0, GL_AMBIENT, lightKa);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightKd);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightKs);
    // position the light
    float lightPos[4] = {light_xpos, light_ypos, light_zpos, 1}; // positional light
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
    glEnable(GL_LIGHT0);                        // MUST enable each light source after configuration
}

///////////////////////////////////////////////////////////////////////////////
// set camera position and lookat direction
///////////////////////////////////////////////////////////////////////////////
void setCamera(float posX, float posY, float posZ, float targetX, float targetY, float targetZ)
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(posX, posY, posZ, targetX, targetY, targetZ, 0, 1, 0); // eye(x,y,z), focal(x,y,z), up(x,y,z)
}

///////////////////////////////////////////////////////////////////////////////
// set the projection matrix as perspective
///////////////////////////////////////////////////////////////////////////////
void toPerspective()
{
    // set viewport to be the entire window
    glViewport(0, 0, (GLsizei)screenWidth, (GLsizei)screenHeight);
    // set perspective viewing frustum
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0f, (float)(screenWidth)/screenHeight, 1.0f, 100.0f); // FOV, AspectRatio, NearClip, FarClip
    // switch to modelview matrix in order to set scene
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

///////////////////////////////////////////////////////////////////////////////
// return true if point is within bounds of room
///////////////////////////////////////////////////////////////////////////////
bool isValidPoint(int x, int y)
{
  if ((x < SCREEN_WIDTH/3 && y < SCREEN_HEIGHT/3) || 
      (x > 2*SCREEN_WIDTH/3 && y < SCREEN_HEIGHT/3) || 
      (x < SCREEN_WIDTH/3 && y > 2*SCREEN_HEIGHT/3) || 
      (x > 2*SCREEN_WIDTH/3 && y > 2*SCREEN_HEIGHT/3)) 
  {
    return false;
  }
  
  return true;
}

///////////////////////////////////////////////////////////////////////////////
// initialize GLUT for windowing
///////////////////////////////////////////////////////////////////////////////
int initGLUT(int argc, char **argv)
{
    // GLUT stuff for windowing
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);   // display mode
    glutInitWindowSize(screenWidth, screenHeight);  // window size
    glutInitWindowPosition(100, 100);               // window location
    // finally, create a window with openGL context
    // Window will not displayed until glutMainLoop() is called
    // it returns a unique ID
    int handle = glutCreateWindow(argv[0]);         // param is the title of window

    // register GLUT callback functions
    glutDisplayFunc(displayCB);
    glutTimerFunc(33, timerCB, 33);                 // redraw only every given millisec
    glutReshapeFunc(reshapeCB);
    glutMouseFunc(mouseCB);
    glutMotionFunc(mouseMotionCB);

    return handle;
}

///////////////////////////////////////////////////////////////////////////////
// initialize GLUT for draw window
///////////////////////////////////////////////////////////////////////////////
int initGLUTDraw(int argc, char **argv)
{
    int drawHandle = glutCreateWindow("Draw Shadow");
    glutSetWindow(drawHandle);
    glutDisplayFunc(displayDrawCB); 
    glutTimerFunc(30, timerDrawCB, 30);
    glutReshapeFunc(reshapeCB);
    glutMouseFunc(mouseDrawCB);
    glutMotionFunc(mouseDrawMotionCB);
    
    return drawHandle;
}
///////////////////////////////////////////////////////////////////////////////
// initialize OpenGL
// disable unused features
///////////////////////////////////////////////////////////////////////////////
void initGL()
{

    glShadeModel(GL_SMOOTH);                    // shading mathod: GL_SMOOTH or GL_FLAT
    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);      // 4-byte pixel alignment
    // enable /disable features
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    //glEnable(GL_CULL_FACE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);


    glClearColor(0, 0, 0, 0);                   // background color
    glClearStencil(0);                          // clear stencil buffer
    glClearDepth(1.0f);                         // 0 is near, 1 is far
    glDepthFunc(GL_LEQUAL);
    initLights();
}

//=============================================================================
// CALLBACKS
//=============================================================================
void timerCB(int millisec)
{
    glutTimerFunc(millisec, timerCB, millisec);
    glutPostRedisplay();
}

void timerDrawCB(int millisec)
{
    glutTimerFunc(millisec, timerDrawCB, millisec);
    // TODO: update lamp info with pixels added/removed
}

void reshapeCB(int w, int h)
{
    screenWidth = w;
    screenHeight = h;
    toPerspective();
}

void mouseCB(int button, int state, int x, int y)
{
    mouseX = x;
    mouseY = y;
    if(button == GLUT_LEFT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseLeftDown = true;
        }
        else if(state == GLUT_UP)
            mouseLeftDown = false;
    }
    else if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseRightDown = true;
        }
        else if(state == GLUT_UP)
            mouseRightDown = false;
    }
}

void mouseMotionCB(int x, int y)
{
    if(mouseLeftDown)
    {
        cameraAngleY += (x - mouseX);
        cameraAngleX += (y - mouseY);
        mouseX = x;
        mouseY = y;
    }
    if(mouseRightDown)
    {
        cameraDistance -= (y - mouseY) * 0.2f;
        mouseY = y;
    }
}

void mouseDrawCB(int button, int state, int x, int y)
{
    drawX = x;
    drawY = y;
    if(button == GLUT_LEFT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            currShadowPixels.clear();
            currShadowPixels.insert(std::make_pair(-1,-1));
            mouseDrawLeftDown = true;
        }
        else if(state == GLUT_UP)
        {
            drawX = x;
            drawY = y;
            for( int i = 0; i < brushWidth; i++){
              for( int j = 0; j < brushWidth; j++){
                if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j))
                    currShadowPixels.insert(std::make_pair(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j));
                    convertedShadowPixels.insert(Vector3f(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j, 5));
              }
            }
            glutPostRedisplay();
            shadowPixels.push_back(currShadowPixels);
            mouseDrawLeftDown = false;
            
            glutSetWindow(mainWindow);
            glutPostRedisplay();
            glutSetWindow(drawWindow);
            //glColor3f(1.0f, 0.0f, 0.0f); 
            //glBegin(GL_POINTS);
            //List::iterator it;
            //Point pixel;
            //for(it = currShadowPixels.begin(); it != currShadowPixels.end(); ++it){
            //  pixel = *it;
            //  glVertex2f(pixel.first,pixel.second);
            //}
            //glEnd();
        }
    }
    else if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            currShadowPixels.clear();
            currShadowPixels.insert(std::make_pair(-2,-2));
            mouseDrawRightDown = true;
        }
        else if(state == GLUT_UP){
            drawX = x;
            drawY = y;
            for( int i = 0; i < brushWidth; i++){
              for( int j = 0; j < brushWidth; j++){
                if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j))
                    currShadowPixels.insert(std::make_pair(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j));
                    convertedShadowPixels.insert(Vector3f(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j, 5));
              }
            }
            glutPostRedisplay();
            mouseDrawRightDown = false;
            shadowPixels.push_back(currShadowPixels);

            //glColor3f(0.0f, 1.0f, 0.0f); 
            //glBegin(GL_POINTS);
            //List::iterator it;
            //Point pixel;
            //for(it = currShadowPixels.begin(); it != currShadowPixels.end(); ++it){
            //  pixel = *it;
            //  glVertex2f(pixel.first,pixel.second);
            //}
            //glEnd();
        }
    }
}

void mouseDrawMotionCB(int x, int y)
{
    if(mouseDrawLeftDown)
    {
        drawX = x;
        drawY = y;
        for( int i = 0; i < brushWidth; i++){
          for( int j = 0; j < brushWidth; j++){
            if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j))
                currShadowPixels.insert(std::make_pair(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j));
                convertedShadowPixels.insert(Vector3f(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j, 5));
          }
        }
        glutPostRedisplay();
    }
    if(mouseDrawRightDown)
    {
        drawX = x;
        drawY = y;
        for( int i = 0; i < brushWidth; i++){
          for( int j = 0; j < brushWidth; j++){
            if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j))
                currShadowPixels.insert(std::make_pair(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j));
          }
        }
        glutPostRedisplay();
    }
}

///////////////////////////////////////////////////////////////////////////////
// Destroy a VBO
// If VBO id is not valid or zero, then OpenGL ignores it silently.
///////////////////////////////////////////////////////////////////////////////
void deleteVBO(const GLuint vboId)
{
    glDeleteBuffers(1, &vboId);
}

///////////////////////////////////////////////////////////////////////////////
// clean up global vars
///////////////////////////////////////////////////////////////////////////////
void clearSharedMem()
{
    deleteVBO(vbo_vertices);
    deleteVBO(vbo_normals);
    deleteVBO(ibo_elements);
}

void exitCB()
{
    clearSharedMem();
}


