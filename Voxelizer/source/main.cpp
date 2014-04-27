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
#include "../include/interface.h"
using namespace std;
using std::stringstream;
using std::cout;
using std::endl;
using std::ends;

#define WALL_THICKNESS		0.05f

/////////////
// GLOBALS //
/////////////

// Voxelized Lamp Mesh Data
GLfloat *lamp_vertices;
GLfloat *lamp_normals;
GLushort *lamp_triangles;

GLuint vbo_vertices;
GLuint vbo_normals;
GLuint ibo_elements;

// Position of the lamp (center) in the space
float lamp_xpos;
float lamp_ypos;

// Position of the lamp (center) in the space
float light_xpos;
float light_ypos;
float light_zpos;

// Room dimension
float room_dim;

unsigned int vertices_size;
unsigned int triangles_size;

// Voxel Lamp Representation for 3d printing 
CompFab::VoxelGrid *g_voxelGrid;

// Used by voxelizer from assn 0;
typedef std::vector<CompFab::Triangle> TriangleList;
TriangleList g_voxelTriangles;
unsigned int voxelRes;

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
    for(unsigned int i = 0; i < g_voxelTriangles.size(); i++){
        CompFab::Triangle triangle = g_voxelTriangles[i];
        if(rayTriangleIntersection(vRay, triangle) == 1){
            numHits ++;
        }
    }
    return numHits;
}

/*
  Load Mesh and construct the voxel grid
  @set g_voxelGrid
*/
bool loadMesh(char *filename, unsigned int dim)
{
    g_voxelTriangles.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1,v2,v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_voxelTriangles.push_back(CompFab::Triangle(v1,v2,v3));
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
  Convert the voxel representation of the lamp into a mesh.
  @sets lamp_vertices, lamp_normals, lamp_triangles
*/
void triangulateVoxelGrid(const char * outfile)
{
    cout << "Trianglulating\n";

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
    // Compute the normals
    mout.compute_norm();
    mout.save_obj(outfile);

    GLfloat p1, p2, p3;
    vertices_size = mout.v.size();
    triangles_size = mout.t.size();

    //Populate the vertices and upload data
    lamp_vertices = new GLfloat[vertices_size*3];
    cout << "parsed triangles size: " << triangles_size << " " << "\n";
    cout << "parsed vertices size: " << vertices_size << " " << "\n";
    for(unsigned int vert =0; vert<mout.v.size(); ++vert)
    {
        p1 = (GLfloat) mout.v[vert][0];
        p2 = (GLfloat) mout.v[vert][1];
        p3 = (GLfloat) mout.v[vert][2];
        lamp_vertices[vert*3] = p1;
        lamp_vertices[vert*3 + 1] = p2;
        lamp_vertices[vert*3 + 2] = p3;
    }

    glGenBuffers(1, &vbo_vertices);                        // create a vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);                    // activate vbo id to use
    glBufferData(GL_ARRAY_BUFFER, vertices_size*3*sizeof(GLfloat), lamp_vertices, GL_STREAM_DRAW); // upload data to video card

    //Populate the normals and upload data
    lamp_normals = new GLfloat[(mout.n.size()*3)];
    for(unsigned int norm =0; norm<mout.n.size(); ++norm)
    {
        p1 = (GLfloat) mout.v[norm][0];
        p2 = (GLfloat) mout.v[norm][1];
        p3 = (GLfloat) mout.v[norm][2];
        lamp_normals[norm*3] = p1;
        lamp_normals[norm*3 + 1] = p2;
        lamp_normals[norm*3 + 2] = p3;
    }
    glGenBuffers(1, &vbo_normals);                        // create a vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);                    // activate vbo id to use
    glBufferData(GL_ARRAY_BUFFER, vertices_size*3*sizeof(GLfloat), lamp_normals, GL_STREAM_DRAW); // upload data to video card

    //Populate the triangle indices and upload data
    lamp_triangles = new GLushort[(mout.t.size()*3)];
    for(unsigned int tri =0; tri<mout.t.size(); ++tri)
    {
        p1 = (GLuint) mout.t[tri][0];
        p2 = (GLuint) mout.t[tri][1];
        p3 = (GLuint) mout.t[tri][2];
        lamp_triangles[tri*3] = p1;
        lamp_triangles[tri*3 + 1] = p2;
        lamp_triangles[tri*3 + 2] = p3;
    }
    glGenBuffers(1, &ibo_elements);                        // create a vbo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);                    // activate vbo id to use
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles_size*3*sizeof(GLushort), lamp_triangles, GL_STREAM_DRAW); // upload data to video card
}

/*
  Voxelize input mesh, update voxel representation, update voxelized mesh representation 
  @sets g_voxelgrid
  @sets lamp_vertices, lamp_normals, lamp_triangles
*/
void voxelizer(char* filename, char* outfilename, unsigned int voxelres) 
{

    unsigned int dim = voxelres; //dimension of voxel grid (e.g. 32x32x32)

    // Construct the voxel grid.
    loadMesh(filename, dim);

    
    // Carve the voxel grid to produce the input lamp voxel representation
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

    // Now produce the lamp mesh representation from the voxel representation
    triangulateVoxelGrid(outfilename);
}

/*
  Create the Scene Data apart from Lamp
  @set light_xpos, light_ypos
TODO:
*/
void createSceneData(float roomDim, float lightXPos, float lightYPos, float lightZPos){ 
    light_xpos = lightXPos;
    light_ypos = lightYPos;
    light_zpos = lightZPos;
    room_dim = roomDim;
}

///////////////
// RENDERING //
///////////////
// REST OF RENDER CODE IS IN interface.cpp AND interface.h //

/*
  Draw room centered at (0,0,0)
*/
void room()
{
	/* ceiling */
	glPushMatrix();
	glTranslatef(0,room_dim/2,0);
	glScalef(room_dim, WALL_THICKNESS, room_dim);
	glutSolidCube( 1.0 );
	glPopMatrix();
	
	/* floor */
	glPushMatrix();
	glTranslatef(0,-room_dim/2,0);
	glScalef(room_dim, WALL_THICKNESS, room_dim);
	glutSolidCube( 1.0 );
	glPopMatrix();
	
	/* right wall */
	glPushMatrix();
	glTranslatef(room_dim/2,0,0);
	glScalef(WALL_THICKNESS, room_dim, room_dim);
	glutSolidCube( 1.0 );
	glPopMatrix();
	
	/* left wall */
	glPushMatrix();
	glTranslatef(-room_dim/2,0,0);
	glScalef(WALL_THICKNESS, room_dim, room_dim);
	glutSolidCube( 1.0 );
	glPopMatrix();
	
	/* back wall */
	glPushMatrix();
	glTranslatef(0, 0, -room_dim/2);
	glScalef(room_dim, room_dim, WALL_THICKNESS);
	glutSolidCube( 1.0 );
	glPopMatrix();
}

/* 
  Callback that actually draws the mesh data 
*/
void displayCB()
{
    // clear buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    setCamera(room_dim*0.5f, room_dim*0.5f, room_dim, room_dim*0.5f, room_dim*0.5f,  room_dim*0.5f);
    // save the initial ModelView matrix before modifying ModelView matrix
    glPushMatrix();
    // tramsform camera
    glTranslatef(0, 0, cameraDistance);
    glRotatef(cameraAngleX, 1, 0, 0);   // pitch
    glRotatef(cameraAngleY, 0, 1, 0);   // heading

    // Enable this for mesh drawing
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // Draw the lamp////////////////////////////////////////////////////////////
    // Set vertex data
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    
    // Set normal data
    glEnableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);
    glNormalPointer(GL_FLOAT, 0, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);
    glDrawElements(GL_TRIANGLES, triangles_size*3, GL_UNSIGNED_SHORT, 0);

    // Draw the cubic room//////////////////////////////////////////////////////
    room();

    glPopMatrix();
    glutSwapBuffers();
}

void displayDrawCB()
{
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho (0, SCREEN_WIDTH, SCREEN_HEIGHT, 0, 0, 1);
	glMatrixMode (GL_MODELVIEW);
	glDisable(GL_DEPTH_TEST);
    //glClear(GL_COLOR_BUFFER_BIT);
    
	glBegin(GL_POINTS);
	glVertex2f(drawX, drawY);
	glEnd();
	
    glutSwapBuffers();
}
//////////
// MAIN //
//////////
int main(int argc, char **argv)
{

    if(argc < 8)
    {
        std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename voxelRes roomDim lightX lightY lightZ\n";
        exit(0);
    }
    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    std::cout<<"Voxel resolution : "<<argv[3]<<"\n";
    std::cout<<"Room dimensions: "<<argv[4]<<"\n";
    std::cout<<"Light position: "<< argv[5] << "," << argv[6] << "," << argv[7] <<"\n";
    initSharedMem();
    // init GLUT and GL
    int mainWindow = initGLUT(argc, argv);
    initGL();

    // init GLUT Draw window
    int drawWindow = initGLUTDraw(argc, argv);
    initGL();
    
    glutSetWindow(mainWindow);

    // register exit callback
    atexit(exitCB);

    int bufferSize;
    voxelRes = atoi(argv[3]);

    // Initialize the lamp mesh
    voxelizer(argv[1], argv[2], voxelRes);
    // Initialize the other scene data
    createSceneData(atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));

    glutMainLoop(); /* Start GLUT event-processing loop */
    return 0;
}

