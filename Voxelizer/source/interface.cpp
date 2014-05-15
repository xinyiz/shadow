#include "../include/interface.h"
using namespace std;
int screenWidth;
int screenHeight;
bool mouseLeftDown;
bool mouseRightDown;
bool mouseDrawLeftDown;
bool mouseDrawRightDown;
int mouseX, mouseY;
int drawX, drawY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance;
int brushWidth;
int mainWindow;
int drawWindow;
int sampleFactor;
int stride;
Queue<List> shadowPixels;
List currShadowPixels;
std::set<Vector3f> pointsToDraw;
Vector3f addIndicator =  Vector3f(-1,-1,-1);
Vector3f remIndicator =  Vector3f(-2,-2,-2);

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
    sampleFactor = 2;
    stride = brushWidth/sampleFactor;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
// initialize lights
///////////////////////////////////////////////////////////////////////////////
void initLights()
{
    // set up light colors (ambient, diffuse, specular)
    GLfloat lightKd[] = {1.0f, 1.0f, 1.0f, 1.0f};  // diffuse light
    GLfloat lightKs[] = {1, 1, 1, 1};           // specular light
    //glLightfv(GL_LIGHT0, GL_AMBIENT, lightKa);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightKd);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightKs);
    // position the light
    float lightPos[4] = {0.5f*room_dim, room_dim, 0.5f*room_dim, 1}; // positional light
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
    glutKeyboardFunc(keyboardCB);

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
    glutKeyboardFunc(keyboardCB);
    
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


    glClearColor(0.7f, 0.7f, 0.7f, 1.0f);                   // background color
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
Vector3f convertTo3DPoint(int x, int y){
      if(x == -1 && y == -1){
        return addIndicator;
      }
      if(x == -2 && y == -2){
        return remIndicator;
      }
      
      // oriented as if the four boxes folded up to become the walls
      if(x >= SCREEN_WIDTH/3 && x <= 2*SCREEN_WIDTH/3 && y >= SCREEN_HEIGHT/3 && y <= 2*SCREEN_HEIGHT/3){
        // middle box = floor of room
        return Vector3f((x-SCREEN_WIDTH/3)*(room_dim/(SCREEN_WIDTH/3)), 0.01, (y-SCREEN_HEIGHT/3)*(room_dim/(SCREEN_HEIGHT/3)));
      } else if (x >= SCREEN_WIDTH/3 && x <= 2*SCREEN_WIDTH/3 && y >=0 && y <= SCREEN_HEIGHT/3) {
        // top box = back wall
       return Vector3f((x-SCREEN_WIDTH/3)*(room_dim/(SCREEN_WIDTH/3)), (SCREEN_HEIGHT/3 - y)*(room_dim/(SCREEN_HEIGHT/3)), 0.03);
      } else if (x >= SCREEN_WIDTH/3 && x <= 2*SCREEN_WIDTH/3 && y >= 2*SCREEN_HEIGHT/3 && y <= SCREEN_HEIGHT){
        // bottom box = front wall
        return Vector3f((x-SCREEN_WIDTH/3)*(room_dim/(SCREEN_WIDTH/3)), (y - SCREEN_HEIGHT)*(room_dim/(SCREEN_HEIGHT/3)) + room_dim, room_dim);
      } else if(x >= 0 && x <= SCREEN_WIDTH/3 && y >= SCREEN_HEIGHT/3 && y <= 2*SCREEN_HEIGHT/3) {
        // left box = left wall
        return Vector3f(0.01,(SCREEN_WIDTH/3 - x)*(room_dim/(SCREEN_WIDTH/3)), (y - SCREEN_HEIGHT/3)*(room_dim/(SCREEN_HEIGHT/3)));
      } else {
        // right box = right wall
        return Vector3f(room_dim-0.1,(x - 2*SCREEN_WIDTH/3)*(room_dim/(SCREEN_WIDTH/3)), (y - SCREEN_HEIGHT/3)*(room_dim/(SCREEN_HEIGHT/3)));
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
            //cout << "Pushing first add point \n";
            currShadowPixels.insert(convertTo3DPoint(-1,-1));
            mouseDrawLeftDown = true;
        }
        else if(state == GLUT_UP)
        {
            drawX = x;
            drawY = y;
            for( int i = 0; i < brushWidth; i++){
              for( int j = 0; j < brushWidth; j++){
                if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j)){
                    Vector3f converted = convertTo3DPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j);
                    if(i%stride == 0){
                        currShadowPixels.insert(converted);
                    }
                    pointsToDraw.insert(converted);
                }
              }
            }
            glutPostRedisplay();
            //cout << "Pushed add updated of size :" << currShadowPixels.size() << "\n";
            shadowPixels.push(currShadowPixels);
            mouseDrawLeftDown = false;
            
            glutSetWindow(mainWindow);
            glutPostRedisplay();
            glutSetWindow(drawWindow);
        }
    }
    else if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            //cout << "Pushing first rem point \n";
            currShadowPixels.clear();
            currShadowPixels.insert(convertTo3DPoint(-2,-2));
            mouseDrawRightDown = true;
        }
        else if(state == GLUT_UP){
            drawX = x;
            drawY = y;
            for( int i = 0; i < brushWidth; i++){
              for( int j = 0; j < brushWidth; j++){
                if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j)){
                    Vector3f converted = convertTo3DPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j);
                    if(i%stride == 0){
                        currShadowPixels.insert(converted);
                    }
                    std::set<Vector3f>::iterator it;
                    for (it = pointsToDraw.begin(); it != pointsToDraw.end();){
                      if (*it == converted){
                        pointsToDraw.erase(it++);
                        break;
                      } else {
                        ++it;
                      }
                    }
                }
              }
            }
            glutPostRedisplay();
            //cout << "Pushed rem updated of size :" << currShadowPixels.size() << "\n";
            mouseDrawRightDown = false;
            shadowPixels.push(currShadowPixels);

            glutSetWindow(mainWindow);
            glutPostRedisplay();
            glutSetWindow(drawWindow);
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
            if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j)){
                Vector3f converted = convertTo3DPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j);
                if(i%stride == 0){
                    currShadowPixels.insert(converted);
                }
                pointsToDraw.insert(converted);
            }
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
            if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j)){
                Vector3f converted = convertTo3DPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j);
                if(i%stride == 0){
                    currShadowPixels.insert(converted);
                }
                std::set<Vector3f>::iterator it;
                for (it = pointsToDraw.begin(); it != pointsToDraw.end();){
                  if (*it == converted){
                    pointsToDraw.erase(it++);
                    break;
                  } else {
                    ++it;
                  }
                }
            }
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


