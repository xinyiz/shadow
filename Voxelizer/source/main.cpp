//Computational Fabrication Assignment #1
// By David Levin 2014
#define GL_GLEXT_PROTOTYPES
#include <GL/glut.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "../include/CompFab.h"
#include "../include/Mesh.h"
#include "../include/vecmath/vecmath.h"
#include "../include/glext.h"
using namespace std;
using std::stringstream;
using std::cout;
using std::endl;
using std::ends;


/////////////
// GLOBALS //
/////////////
typedef std::vector<CompFab::Triangle> TriangleList;
TriangleList voxelLampMeshTriangles;
CompFab::VoxelGrid *g_voxelGrid;
unsigned int voxelRes;
float lightposx = 0.0f;
float lightposy = 0.0f;
float roomDim = 50.0f;
// TODO: VBO for Lamp Mesh


double detMat(double A11, double A12, double A13,
              double A21, double A22, double A23,
              double A31, double A32, double A33)
{
    return A11*(A22*A33 - A23*A32) - A12*(A21*A33 - A23*A31) + A13*(A21*A32 - A22*A31);
}

////////////////
// INITIALIZE //
////////////////

/*
  Ray-Triangle Intersection
  @returns 1 if triangle and ray intersect, 0 otherwise
*/
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    CompFab::Vec3 e1(triangle.m_v2 - triangle.m_v1);
    CompFab::Vec3 e2(triangle.m_v3 - triangle.m_v1);
    CompFab::Vec3 q = ray.m_direction%e2;
    double a = e1*q;
    if(a > -0.000001 and a < 0.000001){
      return 0;
    }

    double f =  1.0/a;
    CompFab::Vec3 s = ray.m_origin - triangle.m_v1;
    double u = f*(s*q); 
    
    if(u < 0.0){
      return 0;
    }
    
    CompFab::Vec3 r = s%e1;
    double v = f*(ray.m_direction*r); 
    if((v < 0.0) or (u + v > 1.0)){
      return 0;
    }
    double t = f*(e2*r);
    if( t <= 0.0 ){
      return 0;
    }
    return 1;

}

/*
  Num Intersections with Ray
  @returns number of intersections with surface made by a ray originating at voxel and cast in direction.
*/
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    
    unsigned int numHits = 0;
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */
    CompFab::RayStruct vRay = CompFab::RayStruct(voxelPos, dir);
    for(unsigned int i = 0; i < voxelLampMeshTriangles.size(); i++){
        CompFab::Triangle triangle = voxelLampMeshTriangles[i];
        if(rayTriangleIntersection(vRay, triangle) == 1){
            numHits ++;
        }
    }
    return numHits;
}

/*
  Load Mesh 
  @set voxelLampMeshTriangles TODO: voxelLampMeshTriangles not voxelized mesh, modify saveVoxelsToObj to do the correct thing 
  @set g_voxelGrid
*/
bool loadMesh(char *filename, unsigned int dim)
{
    voxelLampMeshTriangles.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1,v2,v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        voxelLampMeshTriangles.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    
    return true;
   
}

/*
  Save Voxels to Object File 
  TODO: Modify this to convert voxel representation to triangle rep, store in voxelLampMeshTriangles
*/
void saveVoxelsToObj(const char * outfile)
{
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                  continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}

/*
  Voxelize input mesh, update voxel representation, update voxelized mesh representation 
  TODO: program this. 
*/
void voxelizer(char* filename, unsigned int voxelres) 
{

    unsigned int dim = voxelres; //dimension of voxel grid (e.g. 32x32x32)

    loadMesh(filename, dim);
    

    
    // Cast ray, check if voxel is inside or outside
    // even number of surface intersections = outside (OUT then IN then OUT)
    // odd number = inside (IN then OUT)

    CompFab::Vec3 direction(1.0,0.0,0.0);

    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    CompFab::Vec3 left = g_voxelGrid->m_lowerLeft;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    // Iterate over all voxels in g_voxelGrid and test whether they are inside our outside 
    // of the  surface defined by the triangles in voxelLampMeshTriangles
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                CompFab::Vec3 vPos(left.m_x + ((double)ii)*spacing, left.m_y + ((double)jj)*spacing, left.m_z +((double)kk)*spacing);
                if(numSurfaceIntersections(vPos, direction) % 2 != 0){
                    g_voxelGrid->isInside(ii,jj,kk) = 1;
                }
            }
        }
    }
}


/*
  Create the Scene Data
  @set roomMesh
  @set light
TODO:
*/
void createSceneData(){ 


}
//int main(int argc, char **argv)
//{
//      if(argc < 4)
//      {
//          std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename voxelRes\n";
//          exit(0);
//      }
//      std::cout<<"Load Mesh : "<<argv[1]<<"\n";
//      std::cout<<"Voxel resolution : "<<argv[2]<<"\n";
//
//      voxelRes = atoi(argv[3]); 
//      voxelizer(argv[1], voxelRes);
//
//	    //loadInput();
//	    //glutInit(&argc,argv);
//	    //// We're going to animate it, so double buffer
//	    //glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
//	    //// Initial parameters for window position and size
//	    //glutInitWindowPosition( 60, 60 );
//	    //glutInitWindowSize( 360, 360 );
//	    //glutCreateWindow("Assignment 0");
//	    //// Initialize OpenGL parameters.
//	    //initRendering();
//	    //// Set up callback functions for key presses
//	    //glutKeyboardFunc(keyboardFunc); // Handles "normal" ascii symbols
//	    //glutSpecialFunc(specialFunc);   // Handles "special" keyboard keys
//	    //// Set up the callback function for resizing windows
//	    //glutReshapeFunc( reshapeFunc );
//	    //// Call this whenever window needs redrawing
//	    //glutDisplayFunc( drawScene );
//	    //// Start the main loop.  glutMainLoop never returns.
//	    //glutMainLoop( );
//
//
//	return 0;    // This line is never reached.
//
//      // TODO: for testing right now
//      saveVoxelsToObj(argv[2]);
//      delete g_voxelGrid;
//}



// GLUT CALLBACK functions
void displayCB();
void reshapeCB(int w, int h);
void timerCB(int millisec);
void idleCB();
void keyboardCB(unsigned char key, int x, int y);
void mouseCB(int button, int stat, int x, int y);
void mouseMotionCB(int x, int y);

// CALLBACK function when exit() called ///////////////////////////////////////

void exitCB();
void initGL();
int  initGLUT(int argc, char **argv);
bool initSharedMem();
void clearSharedMem();
void initLights();
void setCamera(float posX, float posY, float posZ, float targetX, float targetY, float targetZ);
GLuint createVBO(const void* data, int dataSize, GLenum target=GL_ARRAY_BUFFER_ARB, GLenum usage=GL_STATIC_DRAW_ARB);
void deleteVBO(const GLuint vboId);
void showInfo();
void toOrtho();
void toPerspective();

// constants
const int   SCREEN_WIDTH    = 400;
const int   SCREEN_HEIGHT   = 300;
const float CAMERA_DISTANCE = 5.0f;
const int   TEXT_WIDTH      = 8;
const int   TEXT_HEIGHT     = 13;

// global variables

void *font = GLUT_BITMAP_8_BY_13;

GLuint vboId = 0;                   // ID of VBO for vertex arrays

int screenWidth;
int screenHeight;
bool mouseLeftDown;
bool mouseRightDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance;
int drawMode = 0;





// cube ///////////////////////////////////////////////////////////////////////
//    v6----- v5
//   /|      /|
//  v1------v0|
//  | |     | |
//  | |v7---|-|v4
//  |/      |/
//  v2------v3

// vertex coords array for glDrawArrays() =====================================
// A cube has 6 sides and each side has 2 triangles, therefore, a cube consists
// of 36 vertices (6 sides * 2 tris * 3 vertices = 36 vertices). And, each
// vertex is 3 components (x,y,z) of floats, therefore, the size of vertex
// array is 108 floats (36 * 3 = 108).
GLfloat vertices[]  = { 1, 1, 1,  -1, 1, 1,  -1,-1, 1,      // v0-v1-v2 (front)
                       -1,-1, 1,   1,-1, 1,   1, 1, 1,      // v2-v3-v0
                        1, 1, 1,   1,-1, 1,   1,-1,-1,      // v0-v3-v4 (right)
                        1,-1,-1,   1, 1,-1,   1, 1, 1,      // v4-v5-v0
                        1, 1, 1,   1, 1,-1,  -1, 1,-1,      // v0-v5-v6 (top)
                       -1, 1,-1,  -1, 1, 1,   1, 1, 1,      // v6-v1-v0
                       -1, 1, 1,  -1, 1,-1,  -1,-1,-1,      // v1-v6-v7 (left)
                       -1,-1,-1,  -1,-1, 1,  -1, 1, 1,      // v7-v2-v1
                       -1,-1,-1,   1,-1,-1,   1,-1, 1,      // v7-v4-v3 (bottom)
                        1,-1, 1,  -1,-1, 1,  -1,-1,-1,      // v3-v2-v7
                        1,-1,-1,  -1,-1,-1,  -1, 1,-1,      // v4-v7-v6 (back)
                       -1, 1,-1,   1, 1,-1,   1,-1,-1 };    // v6-v5-v4
// normal array
GLfloat normals[]   = { 0, 0, 1,   0, 0, 1,   0, 0, 1,      // v0-v1-v2 (front)
                        0, 0, 1,   0, 0, 1,   0, 0, 1,      // v2-v3-v0
                        1, 0, 0,   1, 0, 0,   1, 0, 0,      // v0-v3-v4 (right)
                        1, 0, 0,   1, 0, 0,   1, 0, 0,      // v4-v5-v0
                        0, 1, 0,   0, 1, 0,   0, 1, 0,      // v0-v5-v6 (top)
                        0, 1, 0,   0, 1, 0,   0, 1, 0,      // v6-v1-v0
                       -1, 0, 0,  -1, 0, 0,  -1, 0, 0,      // v1-v6-v7 (left)
                       -1, 0, 0,  -1, 0, 0,  -1, 0, 0,      // v7-v2-v1
                        0,-1, 0,   0,-1, 0,   0,-1, 0,      // v7-v4-v3 (bottom)
                        0,-1, 0,   0,-1, 0,   0,-1, 0,      // v3-v2-v7
                        0, 0,-1,   0, 0,-1,   0, 0,-1,      // v4-v7-v6 (back)
                        0, 0,-1,   0, 0,-1,   0, 0,-1 };    // v6-v5-v4
// color array
GLfloat colors[]    = { 1, 1, 1,   1, 1, 0,   1, 0, 0,      // v0-v1-v2 (front)
                        1, 0, 0,   1, 0, 1,   1, 1, 1,      // v2-v3-v0
                        1, 1, 1,   1, 0, 1,   0, 0, 1,      // v0-v3-v4 (right)
                        0, 0, 1,   0, 1, 1,   1, 1, 1,      // v4-v5-v0
                        1, 1, 1,   0, 1, 1,   0, 1, 0,      // v0-v5-v6 (top)
                        0, 1, 0,   1, 1, 0,   1, 1, 1,      // v6-v1-v0
                        1, 1, 0,   0, 1, 0,   0, 0, 0,      // v1-v6-v7 (left)
                        0, 0, 0,   1, 0, 0,   1, 1, 0,      // v7-v2-v1
                        0, 0, 0,   0, 0, 1,   1, 0, 1,      // v7-v4-v3 (bottom)
                        1, 0, 1,   1, 0, 0,   0, 0, 0,      // v3-v2-v7
                        0, 0, 1,   0, 0, 0,   0, 1, 0,      // v4-v7-v6 (back)
                        0, 1, 0,   0, 1, 1,   0, 0, 1 };    // v6-v5-v4

// function pointers for VBO Extension
// Windows needs to get function pointers from ICD OpenGL drivers,
// because opengl32.dll does not support extensions higher than v1.1.

#ifdef _WIN32
PFNGLGENBUFFERSARBPROC            pglGenBuffersARB = 0;             // VBO Name Generation Procedure
PFNGLBINDBUFFERARBPROC            pglBindBufferARB = 0;             // VBO Bind Procedure
PFNGLBUFFERDATAARBPROC            pglBufferDataARB = 0;             // VBO Data Loading Procedure
PFNGLBUFFERSUBDATAARBPROC         pglBufferSubDataARB = 0;          // VBO Sub Data Loading Procedure
PFNGLDELETEBUFFERSARBPROC         pglDeleteBuffersARB = 0;          // VBO Deletion Procedure
PFNGLGETBUFFERPARAMETERIVARBPROC  pglGetBufferParameterivARB = 0;   // return various parameters of VBO
PFNGLMAPBUFFERARBPROC             pglMapBufferARB = 0;              // map VBO procedure
PFNGLUNMAPBUFFERARBPROC           pglUnmapBufferARB = 0;            // unmap VBO procedure

#define glGenBuffersARB           pglGenBuffersARB
#define glBindBufferARB           pglBindBufferARB
#define glBufferDataARB           pglBufferDataARB
#define glBufferSubDataARB        pglBufferSubDataARB
#define glDeleteBuffersARB        pglDeleteBuffersARB
#define glGetBufferParameterivARB pglGetBufferParameterivARB
#define glMapBufferARB            pglMapBufferARB
#define glUnmapBufferARB          pglUnmapBufferARB
#endif







///////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)

{

    initSharedMem();
    // init GLUT and GL
    initGLUT(argc, argv);
    initGL();
    // register exit callback
    atexit(exitCB);
#ifdef _WIN32
    // get pointers to GL functions
    glGenBuffersARB = (PFNGLGENBUFFERSARBPROC)wglGetProcAddress("glGenBuffersARB");
    glBindBufferARB = (PFNGLBINDBUFFERARBPROC)wglGetProcAddress("glBindBufferARB");
    glBufferDataARB = (PFNGLBUFFERDATAARBPROC)wglGetProcAddress("glBufferDataARB");
    glBufferSubDataARB = (PFNGLBUFFERSUBDATAARBPROC)wglGetProcAddress("glBufferSubDataARB");
    glDeleteBuffersARB = (PFNGLDELETEBUFFERSARBPROC)wglGetProcAddress("glDeleteBuffersARB");
    glGetBufferParameterivARB = (PFNGLGETBUFFERPARAMETERIVARBPROC)wglGetProcAddress("glGetBufferParameterivARB");
    glMapBufferARB = (PFNGLMAPBUFFERARBPROC)wglGetProcAddress("glMapBufferARB");
    glUnmapBufferARB = (PFNGLUNMAPBUFFERARBPROC)wglGetProcAddress("glUnmapBufferARB");
#else 
#endif


    int bufferSize;
    // create vertex buffer objects, you need to delete them when program exits
    // Try to put both vertex coords array, vertex normal array and vertex color in the same buffer object.
    // glBufferDataARB with NULL pointer reserves only memory space.
    // Copy actual data with 2 calls of glBufferSubDataARB, one for vertex coords and one for normals.
    // target flag is GL_ARRAY_BUFFER_ARB, and usage flag is GL_STATIC_DRAW_ARB

    glGenBuffersARB(1, &vboId);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboId);
    glBufferDataARB(GL_ARRAY_BUFFER_ARB, sizeof(vertices)+sizeof(normals)+sizeof(colors), 0, GL_STATIC_DRAW_ARB);
    glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, sizeof(vertices), vertices);                             // copy vertices starting from 0 offest
    glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, sizeof(vertices), sizeof(normals), normals);                // copy normals after vertices
    glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, sizeof(vertices)+sizeof(normals), sizeof(colors), colors);  // copy colours after normals

    // the last GLUT call (LOOP)
    // window will be shown and display callback is triggered by events
    // NOTE: this call never return main().
    glutMainLoop(); /* Start GLUT event-processing loop */
    return 0;
}



///////////////////////////////////////////////////////////////////////////////

// initialize GLUT for windowing

///////////////////////////////////////////////////////////////////////////////

int initGLUT(int argc, char **argv)

{

    // GLUT stuff for windowing
    // initialization openGL window.
    // it is called before any other GLUT routine
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
    //glutIdleFunc(idleCB);                           // redraw when idle
    glutReshapeFunc(reshapeCB);
    glutKeyboardFunc(keyboardCB);
    glutMouseFunc(mouseCB);
    glutMotionFunc(mouseMotionCB);
    return handle;
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
    //glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    //glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_CULL_FACE);
     // track material ambient and diffuse from surface color, call it before glEnable(GL_COLOR_MATERIAL)
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glClearColor(0, 0, 0, 0);                   // background color
    glClearStencil(0);                          // clear stencil buffer
    glClearDepth(1.0f);                         // 0 is near, 1 is far
    glDepthFunc(GL_LEQUAL);
    initLights();
}

///////////////////////////////////////////////////////////////////////////////

// initialize global variables

///////////////////////////////////////////////////////////////////////////////

bool initSharedMem()
{
    screenWidth = SCREEN_WIDTH;
    screenHeight = SCREEN_HEIGHT;
    mouseLeftDown = mouseRightDown = false;
    mouseX = mouseY = 0;
    cameraAngleX = cameraAngleY = 0.0f;
    cameraDistance = CAMERA_DISTANCE;
    drawMode = 0; // 0:fill, 1: wireframe, 2:points
    return true;
}
///////////////////////////////////////////////////////////////////////////////

// clean up global vars

///////////////////////////////////////////////////////////////////////////////

void clearSharedMem()
{
    deleteVBO(vboId);

    vboId = 0;
}

///////////////////////////////////////////////////////////////////////////////

// initialize lights

///////////////////////////////////////////////////////////////////////////////

void initLights()

{
    // set up light colors (ambient, diffuse, specular)
    GLfloat lightKa[] = {.2f, .2f, .2f, 1.0f};  // ambient light
    GLfloat lightKd[] = {.7f, .7f, .7f, 1.0f};  // diffuse light
    GLfloat lightKs[] = {1, 1, 1, 1};           // specular light
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightKa);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightKd);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightKs);
    // position the light
    float lightPos[4] = {0, 0, 20, 1}; // positional light
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
// generate vertex buffer object and bind it with its data
// You must give 2 hints about data usage; target and mode, so that OpenGL can
// decide which data should be stored and its location.
// VBO works with 2 different targets; GL_ARRAY_BUFFER_ARB for vertex arrays
// and GL_ELEMENT_ARRAY_BUFFER_ARB for index array in glDrawElements().
// The default target is GL_ARRAY_BUFFER_ARB.
// By default, usage mode is set as GL_STATIC_DRAW_ARB.
// Other usages are GL_STREAM_DRAW_ARB, GL_STREAM_READ_ARB, GL_STREAM_COPY_ARB,
// GL_STATIC_DRAW_ARB, GL_STATIC_READ_ARB, GL_STATIC_COPY_ARB,
// GL_DYNAMIC_DRAW_ARB, GL_DYNAMIC_READ_ARB, GL_DYNAMIC_COPY_ARB.
///////////////////////////////////////////////////////////////////////////////

GLuint createVBO(const void* data, int dataSize, GLenum target, GLenum usage)

{
    GLuint id = 0;  // 0 is reserved, glGenBuffersARB() will return non-zero id if success
    glGenBuffersARB(1, &id);                        // create a vbo
    glBindBufferARB(target, id);                    // activate vbo id to use
    glBufferDataARB(target, dataSize, data, usage); // upload data to video card
    // check data size in VBO is same as input array, if not return 0 and delete VBO
    int bufferSize = 0;
    glGetBufferParameterivARB(target, GL_BUFFER_SIZE_ARB, &bufferSize);
    if(dataSize != bufferSize)
    {
        glDeleteBuffersARB(1, &id);
        id = 0;
        std::cout << "[createVBO()] Data size is mismatch with input array\n";
    }
    return id;      // return VBO id
}
///////////////////////////////////////////////////////////////////////////////

// destroy a VBO

// If VBO id is not valid or zero, then OpenGL ignores it silently.

///////////////////////////////////////////////////////////////////////////////
void deleteVBO(const GLuint vboId)
{
    glDeleteBuffersARB(1, &vboId);
}

///////////////////////////////////////////////////////////////////////////////

// set projection matrix as orthogonal

///////////////////////////////////////////////////////////////////////////////

void toOrtho()

{
    // set viewport to be the entire window
    glViewport(0, 0, (GLsizei)screenWidth, (GLsizei)screenHeight);
    // set orthographic viewing frustum
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, screenWidth, 0, screenHeight, -1, 1);
    // switch to modelview matrix in order to set scene
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
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
    gluPerspective(60.0f, (float)(screenWidth)/screenHeight, 1.0f, 1000.0f); // FOV, AspectRatio, NearClip, FarClip
    // switch to modelview matrix in order to set scene
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}




//=============================================================================

// CALLBACKS

//=============================================================================



void displayCB()

{
    // clear buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    // save the initial ModelView matrix before modifying ModelView matrix
    glPushMatrix();
    // tramsform camera
    glTranslatef(0, 0, -cameraDistance);
    glRotatef(cameraAngleX, 1, 0, 0);   // pitch
    glRotatef(cameraAngleY, 0, 1, 0);   // heading
    // bind VBOs with IDs and set the buffer offsets of the bound VBOs
    // When buffer object is bound with its ID, all pointers in gl*Pointer()
    // are treated as offset instead of real pointer.
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, vboId);
    // enable vertex arrays
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);
    // before draw, specify vertex and index arrays with their offsets
    glNormalPointer(GL_FLOAT, 0, (void*)sizeof(vertices));
    glColorPointer(3, GL_FLOAT, 0, (void*)(sizeof(vertices)+sizeof(normals)));
    glVertexPointer(3, GL_FLOAT, 0, 0);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glDisableClientState(GL_VERTEX_ARRAY);  // disable vertex arrays
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    // it is good idea to release VBOs with ID 0 after use.
    // Once bound with 0, all pointers in gl*Pointer() behave as real
    // pointer, so, normal vertex array operations are re-activated
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
    glPopMatrix();
    glutSwapBuffers();
}

void reshapeCB(int w, int h)
{
    screenWidth = w;
    screenHeight = h;
    toPerspective();
}

void timerCB(int millisec)
{
    glutTimerFunc(millisec, timerCB, millisec);
    glutPostRedisplay();
}

void idleCB()
{
    glutPostRedisplay();
}

void keyboardCB(unsigned char key, int x, int y)
{
    switch(key)
    {
    case 27: // ESCAPE
        exit(0);
        break;
    case 'd': // switch rendering modes (fill -> wire -> point)
    case 'D':
        drawMode = ++drawMode % 3;
        if(drawMode == 0)        // fill mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_CULL_FACE);
        }
        else if(drawMode == 1)  // wireframe mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_CULL_FACE);
        }
        else                    // point mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_CULL_FACE);
        }
        break;
    default:
        ;
    }
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

void exitCB()
{
    clearSharedMem();
}
