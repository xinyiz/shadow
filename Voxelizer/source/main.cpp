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
#include <assert.h>
#include <ctime>
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

///////////////
// CONSTANTS //
///////////////
int c_res = 10;
int c_stride = 360/c_res;
int nextVoxelLookup[36] = 
{ 
     0, 0,-1, 
     0, 0,-1, 
     0,-1, 0, 
     0,-1, 0, 
     1, 0, 0, 
     1, 0, 0, 
     0, 1, 0, 
     0, 1, 0, 
    -1, 0, 0, 
    -1, 0, 0, 
     0, 0, 1, 
     0, 0, 1, 
};


inline void removeActiveTriangles(unsigned int start);
inline void addActiveTriangles(unsigned int start);
inline void updateNextVoxel(unsigned int &curr_i, unsigned int &curr_j, unsigned int &curr_k, unsigned int tri);
void voxelsIntersect(int ii, int jj, int kk, CompFab::Vec3 &shadePoint, bool add);
/////////////
// GLOBALS //
/////////////

typedef std::vector<CompFab::Triangle> TriangleList;

/* Voxelized Lamp Mesh Data for 3D PRINTING */
// State of carved lamp represented by g_carvedLampMesh->activeTriangles, 
// indices of the untouched triangles in the lamp. 
//TODO: update saved obj file to only include active triangles
Mesh g_carvedLampMesh;

/* Voxelized Lamp Mesh Data for RENDERING */
GLfloat *lamp_vertices;
GLfloat *lamp_normals;
GLushort *lamp_triangles;
GLushort *lamp_triangles_test;
int numAct = 0;
int numInact = 0;

GLuint vbo_vertices;
GLuint vbo_normals;
GLuint ibo_elements;
GLuint ibo_elements_test;

/* Voxel Data Structure for Updating Lamp */
CompFab::VoxelGrid *g_lampVoxelGrid;
double gridSpacing;
CompFab::Vec3 gridLLeft;
int nx, ny, nz;

// Triangles in the input lamp mesh - used for voxel intersection checking
TriangleList g_inputLampTriangles;

// Used by voxelizer from assn 0;
TriangleList g_voxelTriangles;
unsigned int voxelRes;

/* Scene Data */
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


bool cleared = false;

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
float rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
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
    //cout << "T:" << t << "\n";
    if( t <= 0.0 ){
      return 0;
    }
    return (float)t;
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
        if(rayTriangleIntersection(vRay, triangle)){
            numHits ++;
        }
    }
    return numHits;
}

/*
  Load Mesh and construct the voxel grid
  @set g_lampVoxelGrid
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
    
    g_lampVoxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);
    gridSpacing = g_lampVoxelGrid->m_spacing;
    gridLLeft = g_lampVoxelGrid->m_lowerLeft;
    nx = g_lampVoxelGrid->m_dimX;
    ny = g_lampVoxelGrid->m_dimY;
    nz = g_lampVoxelGrid->m_dimZ;


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
    CompFab::Vec3 hspacing(0.5*gridSpacing, 0.5*gridSpacing, 0.5*gridSpacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_lampVoxelGrid->isInside(ii,jj,kk)){
                  continue;
                }
                CompFab::Vec3 coord(((double)ii)*gridSpacing, ((double)jj)*gridSpacing, ((double)kk)*gridSpacing);
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
    cout << "Triangulating\n";

    Mesh box;
    CompFab::Vec3 hspacing(0.5*gridSpacing, 0.5*gridSpacing, 0.5*gridSpacing);

    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_lampVoxelGrid->isInside(ii,jj,kk)){
                  continue;
                }
                CompFab::Vec3 coord(((double)ii)*gridSpacing, ((double)jj)*gridSpacing, ((double)kk)*gridSpacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                int triIndex = g_carvedLampMesh.append(box);
                //Denote the start index of the 12 triangles for this voxel
                //the indices are contiguous
                if( ii == 10 && jj == 10 && kk == 10){
                  cout << "TRIANGLES FOR BOX 10,10,10: \n";
                  CompFab::Vec3 v1, v2, v3;
                  for(int i = 0; i < 12; i++){
                    v1 = box.v[box.t[i][0]];
                    v2 = box.v[box.t[i][1]];
                    v3 = box.v[box.t[i][2]];
                    cout << "v1: " << v1.m_x << "," << v1.m_y << "," << v1.m_z << "\n";
                    cout << "v2: " << v2.m_x << "," << v2.m_y << "," << v2.m_z << "\n";
                    cout << "v3: " << v3.m_x << "," << v3.m_y << "," << v3.m_z << "\n";
                    cout << "\n";

                  }
                    cout << "TARGET 1: " << triIndex;
                }
                g_lampVoxelGrid->setTrianglesIndex(ii,jj,kk,triIndex);
            }
        }
    }
    // Compute the normals
    cout << "Computing normals...\n";
    g_carvedLampMesh.compute_norm();
    // Initialize the list of active triangles to be all
    cout << "Initializing active triangles...\n";
    g_carvedLampMesh.init_active_triangles();

    CompFab::Vec3 v1,v2,v3;

    cout << "Inserting into inputLampTriangles...\n";
    for(unsigned int tri =0; tri<g_carvedLampMesh.t.size(); ++tri)
    {
        v1 = g_carvedLampMesh.v[g_carvedLampMesh.t[tri][0]];
        v2 = g_carvedLampMesh.v[g_carvedLampMesh.t[tri][1]];
        v3 = g_carvedLampMesh.v[g_carvedLampMesh.t[tri][2]];
        g_inputLampTriangles.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //cout << "Testing...\n";
    //TODO: testing - remove
    //// jitter the ray
    //srand(time(NULL));
    //float rand1 = static_cast <float> (rand()/ static_cast<float> (RAND_MAX))*0.01f;
    //float rand2 = static_cast <float> (rand()/ static_cast<float> (RAND_MAX))*0.01f;
    //float rand3 = static_cast <float> (rand()/ static_cast<float> (RAND_MAX))*0.01f;
    //cout << "Random nums:" << rand1 << "," << rand2 << "," << rand3;
    //CompFab::Vec3 spoint(5.0+rand1,5.0+rand2,5.0+rand3);
    ////CompFab::Vec3 spoint(5.1,5.3,5.4);
    ////CompFab::Vec3 spoint(gridLLeft.m_x + ((double)3)*gridSpacing + rand1, gridLLeft.m_y + ((double)3)*gridSpacing+rand2, 6 + gridLLeft.m_z +((double)3)*gridSpacing+rand3);
    //voxelsIntersect(voxelRes/2, voxelRes/2, voxelRes/2, spoint, false);

    cout << "Saving...\n";
    g_carvedLampMesh.save_obj(outfile);

    GLfloat p1, p2, p3;
    vertices_size = g_carvedLampMesh.v.size();
    triangles_size = g_carvedLampMesh.t.size();

    //Populate the vertices and upload data
    lamp_vertices = new GLfloat[vertices_size*3];
    for(unsigned int vert =0; vert<g_carvedLampMesh.v.size(); ++vert)
    {
        p1 = (GLfloat) g_carvedLampMesh.v[vert][0];
        p2 = (GLfloat) g_carvedLampMesh.v[vert][1];
        p3 = (GLfloat) g_carvedLampMesh.v[vert][2];
        lamp_vertices[vert*3] = p1;
        lamp_vertices[vert*3 + 1] = p2;
        lamp_vertices[vert*3 + 2] = p3;
    }

    glGenBuffers(1, &vbo_vertices);                        // create a vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);                    // activate vbo id to use
    glBufferData(GL_ARRAY_BUFFER, vertices_size*3*sizeof(GLfloat), lamp_vertices, GL_STREAM_DRAW); // upload data to video card

    //Populate the normals and upload data
    lamp_normals = new GLfloat[(g_carvedLampMesh.n.size()*3)];
    for(unsigned int norm =0; norm<g_carvedLampMesh.n.size(); ++norm)
    {
        p1 = (GLfloat) g_carvedLampMesh.v[norm][0];
        p2 = (GLfloat) g_carvedLampMesh.v[norm][1];
        p3 = (GLfloat) g_carvedLampMesh.v[norm][2];
        lamp_normals[norm*3] = p1;
        lamp_normals[norm*3 + 1] = p2;
        lamp_normals[norm*3 + 2] = p3;
    }
    glGenBuffers(1, &vbo_normals);                        // create a vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);                    // activate vbo id to use
    glBufferData(GL_ARRAY_BUFFER, vertices_size*3*sizeof(GLfloat), lamp_normals, GL_STREAM_DRAW); // upload data to video card

    //Populate the triangle indices and upload data
    lamp_triangles = new GLushort[(g_carvedLampMesh.t.size()*3)];
    lamp_triangles_test = new GLushort[(g_carvedLampMesh.t.size()*3)];
    for(unsigned int tri =0; tri<g_carvedLampMesh.t.size(); ++tri)
    {
        p1 = (GLuint) g_carvedLampMesh.t[tri][0];
        p2 = (GLuint) g_carvedLampMesh.t[tri][1];
        p3 = (GLuint) g_carvedLampMesh.t[tri][2];
      if(g_carvedLampMesh.activeTriangles.find(tri) != g_carvedLampMesh.activeTriangles.end())
      {
        lamp_triangles[numAct*3] = p1;
        lamp_triangles[numAct*3 + 1] = p2;
        lamp_triangles[numAct*3 + 2] = p3;
        numAct+=1;
      } else {
        lamp_triangles_test[numInact*3] = p1;
        lamp_triangles_test[numInact*3 + 1] = p2;
        lamp_triangles_test[numInact*3 + 2] = p3;
        numInact+=1;
      }

    }
    assert((numAct + numInact) == triangles_size);
    glGenBuffers(1, &ibo_elements);                        // create a vbo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);                    // activate vbo id to use
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, numAct*3*sizeof(GLushort), lamp_triangles, GL_STREAM_DRAW); // upload data to video card

    glGenBuffers(1, &ibo_elements_test);                        // create a vbo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements_test);                    // activate vbo id to use
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, numInact*3*sizeof(GLushort), lamp_triangles_test, GL_STREAM_DRAW); // upload data to video card
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

    cout << "m_lowerleft" << gridLLeft.m_x << "," << gridLLeft.m_y << "," << gridLLeft.m_z;
    
    CompFab::Vec3 hspacing(0.5*gridSpacing, 0.5*gridSpacing, 0.5*gridSpacing);
    
    // Iterate over all voxels in g_lampVoxelGrid and test whether they are inside our outside 
    // of the  surface defined by the triangles in voxelLampMeshTriangles
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                CompFab::Vec3 vPos(gridLLeft.m_x + ((double)ii)*gridSpacing, gridLLeft.m_y + ((double)jj)*gridSpacing, gridLLeft.m_z +((double)kk)*gridSpacing);
                if(numSurfaceIntersections(vPos, direction) % 2 != 0){
                    g_lampVoxelGrid->isInside(ii,jj,kk) = 1;
                    g_lampVoxelGrid->isCarved(ii,jj,kk) = 0;
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
//////////////
// UPDATING //
//////////////
/* 
  Maps triangle to the voxel bordering face containing 
  triangle - mapping is change to curr i,j,k voxel index.
*/
//TODO: doesn't consider edge cases


inline void updateNextVoxel(unsigned int &curr_i, unsigned int &curr_j, unsigned int &curr_k, unsigned int tri){
    curr_i = curr_i + nextVoxelLookup[tri*3];
    curr_j = curr_j + nextVoxelLookup[tri*3 + 1];
    curr_k = curr_k + nextVoxelLookup[tri*3 + 2];
}

inline void addActiveTriangles(unsigned int start){
    for(int i = start; i < start+12; i++){
      g_carvedLampMesh.activeTriangles.insert(i);
    }
}

inline void removeActiveTriangles(unsigned int start){
    for(int i = start; i < start+12; i++){
      g_carvedLampMesh.activeTriangles.erase(i);
    }
}
void printTriangle(CompFab::Triangle T){
  cout << "Triangle---\n";
  cout << "v1: " << T.m_v1.m_x << "," << T.m_v1.m_y << "," << T.m_v1.m_z << "\n";
  cout << "v2: " << T.m_v2.m_x << "," << T.m_v2.m_y << "," << T.m_v2.m_z << "\n";
  cout << "v3: " << T.m_v3.m_x << "," << T.m_v3.m_y << "," << T.m_v3.m_z << "\n";
}
/*
  @param ii, jj, kk - voxel containing the light source  
  @param shadePoint - the end point of the light ray
  Correctly updates the lamp mesh information for 3d printing and rendering
*/
void voxelsIntersect(int ii, int jj, int kk, CompFab::Vec3 &shadePoint, bool add){
    cout << "Intersecting ray with lamp for voxel updates...";
    std::vector<int> voxelIndices;
    CompFab::Vec3 vPos(gridLLeft.m_x + ((double)ii)*gridSpacing + gridSpacing*0.5f, gridLLeft.m_y + ((double)jj)*gridSpacing + gridSpacing*0.5f, gridLLeft.m_z +((double)kk)*gridSpacing+ gridSpacing*0.5f);
    CompFab::Vec3 dir = (shadePoint - vPos);
    //CompFab::Vec3 dir(0,-1,0);
    dir.normalize();
    CompFab::RayStruct vRay = CompFab::RayStruct(vPos, dir);
    //cout << "Origin:" <<  vPos.m_x << "," << vPos.m_y << "," << vPos.m_z << "\n";

    //cout << "Shade point:" <<  shadePoint.m_x << "," << shadePoint.m_y << "," << shadePoint.m_z << "\n";
    //cout << "Ray dir:" <<  dir.m_x << "," << dir.m_y << "," << dir.m_z << "\n";

    unsigned int curr_i = ii;
    unsigned int curr_j = jj;
    unsigned int curr_k = kk;
    unsigned int startTriangle;
    float prev_d, curr_d;       // rayTriangleIntersection dist to track entering/exit face


    startTriangle = g_lampVoxelGrid->getFirstTriangle(ii, jj, kk);
    assert(g_lampVoxelGrid->isInside(ii,jj,kk));
    cout << "Start from first voxel..." << ii << "," << jj << "," << kk << " triangle num " << startTriangle << "\n";
    for(unsigned int tri = startTriangle; tri< startTriangle + 12; ++tri)
    {
        printTriangle(g_inputLampTriangles[tri]);
        prev_d = rayTriangleIntersection(vRay, g_inputLampTriangles[tri]);
        //cout << "prev_d... " << prev_d << "\n";
        if(prev_d){
            if(add){
                if(g_lampVoxelGrid->isInside(curr_i, curr_j, curr_k) == 1 && 
                   g_lampVoxelGrid->isCarved(curr_i, curr_j, curr_k) == 1){
                    addActiveTriangles(startTriangle);          //Update g_carvedLampMesh
                    g_lampVoxelGrid->isCarved(curr_i, curr_j, curr_k) = 0;
                    cout << "Adding Triangles start at: " << startTriangle << "\n";
                }
            } else {
                if(g_lampVoxelGrid->isInside(curr_i, curr_j, curr_k) == 1 && 
                   g_lampVoxelGrid->isCarved(curr_i, curr_j, curr_k) == 0){
                    removeActiveTriangles(startTriangle);       //Update g_carvedLampMesh 
                    g_lampVoxelGrid->isCarved(curr_i, curr_j, curr_k) = 1;
                    cout << "Removing Triangles start at: " << startTriangle << "\n";
                }
            }
            updateNextVoxel(curr_i,curr_j,curr_k,(tri % 12));
            break;
        }
    }

    // Case 1:  ADD 
    // Iterate until ray exits lamp bounds 
    if(add){
        cout << "Adding back voxels in ray path...\n";
        while(true){
            if(g_lampVoxelGrid->isInside(curr_i, curr_j, curr_k) == 0)
                return;

            startTriangle = g_lampVoxelGrid->getFirstTriangle(curr_i, curr_j, curr_k);
            for(unsigned int tri = startTriangle; tri< startTriangle + 12; ++tri)
            {
                // Figure out which voxel to examine next
                curr_d = rayTriangleIntersection(vRay, g_inputLampTriangles[tri]);
                if(prev_d < curr_d){    // Check triangle exit face triangle
                    if(g_lampVoxelGrid->isInside(curr_i, curr_j, curr_k) == 1 && 
                       g_lampVoxelGrid->isCarved(curr_i, curr_j, curr_k) == 1){
                        addActiveTriangles(startTriangle);       //Update g_carvedLampMesh 
                        g_lampVoxelGrid->isCarved(curr_i, curr_j, curr_k) = 0;
                        cout << "Adding Triangles start at: " << startTriangle << "\n";
                    }
                    updateNextVoxel(curr_i,curr_j,curr_k,(tri % 12));
                    break;
                }
            }
            prev_d = curr_d;
        }
    } else{ // Case 2:  REMOVE 
        cout << "Removing voxels in ray path...\n";
        while(true){
            cout << "Curr voxel: " << curr_i << "," << curr_j << "," << curr_k << "\n";
            if(g_lampVoxelGrid->isInside(curr_i, curr_j, curr_k) == 0)
                return;

            startTriangle = g_lampVoxelGrid->getFirstTriangle(curr_i, curr_j, curr_k);
            for(unsigned int tri = startTriangle; tri< startTriangle + 12; tri++)
            {
                // Figure out which voxel to examine next

                curr_d = rayTriangleIntersection(vRay, g_inputLampTriangles[tri]);
                //cout << "prev_d curr_d:" <<  prev_d << "," << curr_d << "\n";
                //cout << "start triangle:" << tri << "\n";
                if(prev_d < curr_d){    // Check triangle exit face triangle
                    if(g_lampVoxelGrid->isInside(curr_i, curr_j, curr_k) == 1 && 
                       g_lampVoxelGrid->isCarved(curr_i, curr_j, curr_k) == 0){
                        removeActiveTriangles(startTriangle);       //Update g_carvedLampMesh 
                        g_lampVoxelGrid->isCarved(curr_i, curr_j, curr_k) = 1;
                        cout << "Removing Triangles start at: " << startTriangle << "\n";
                    }
                    updateNextVoxel(curr_i,curr_j,curr_k,(tri % 12));
                    break;
                }
            }
            prev_d = curr_d;
        }
    }
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
    glTranslatef(room_dim/2,room_dim,room_dim/2);
    glScalef(room_dim, WALL_THICKNESS, room_dim);
    glutSolidCube( 1.0 );
    glPopMatrix();
    
    /* floor */
    glPushMatrix();
    glTranslatef(room_dim/2,0,room_dim/2);
    glScalef(room_dim, WALL_THICKNESS, room_dim);
    glutSolidCube( 1.0 );
    glPopMatrix();
	
    /* right wall */
    glPushMatrix();
    glTranslatef(room_dim,room_dim/2,room_dim/2);
    glScalef(WALL_THICKNESS, room_dim, room_dim);
    glutSolidCube( 1.0 );
    glPopMatrix();
	
    /* left wall */
    glPushMatrix();
    glTranslatef(0,room_dim/2,room_dim/2);
    glScalef(WALL_THICKNESS, room_dim, room_dim);
    glutSolidCube( 1.0 );
    glPopMatrix();
   
    /* back wall */
    glPushMatrix();
    glTranslatef(room_dim/2,room_dim/2,0);
    glScalef(room_dim, room_dim, WALL_THICKNESS);
    glutSolidCube( 1.0 );
    glPopMatrix();

}

void updateTriangleVBOs(){
    //Populate the triangle indices and upload data
    GLfloat p1, p2, p3;
    numAct = 0;
    numInact = 0;
    for(unsigned int tri =0; tri<g_carvedLampMesh.t.size(); ++tri)
    {
        p1 = (GLuint) g_carvedLampMesh.t[tri][0];
        p2 = (GLuint) g_carvedLampMesh.t[tri][1];
        p3 = (GLuint) g_carvedLampMesh.t[tri][2];
      if(g_carvedLampMesh.activeTriangles.find(tri) != g_carvedLampMesh.activeTriangles.end())
      {
        lamp_triangles[numAct*3] = p1;
        lamp_triangles[numAct*3 + 1] = p2;
        lamp_triangles[numAct*3 + 2] = p3;
        numAct+=1;
      } else {
        lamp_triangles_test[numInact*3] = p1;
        lamp_triangles_test[numInact*3 + 1] = p2;
        lamp_triangles_test[numInact*3 + 2] = p3;
        numInact+=1;
      }

    }
    assert((numAct + numInact) == triangles_size);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);                    // activate vbo id to use
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, numAct*3*sizeof(GLushort), lamp_triangles, GL_STREAM_DRAW); // upload data to video card

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements_test);                    // activate vbo id to use
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, numInact*3*sizeof(GLushort), lamp_triangles_test, GL_STREAM_DRAW); // upload data to video card

}
void updateLamp(std::set<Vector3f> points){
   cout << "Points size: " << points.size() << "\n";
   bool add = false;

   std::set<Vector3f>::iterator it;
   Vector3f point;
   // first point in set is indicator
   for(it = points.begin(); it != points.end(); ++it){
     point = *it;
     // since points.find() wasn't finding indicators
     if (point == addIndicator) {
        cout << "adding points" << "\n";
        add = true;
     } else if (point == remIndicator) {
        cout << "removing points" << "\n";
        add = false;
     } else {
     // Jitter the shadow point
     srand(time(NULL));
     float rand1 = static_cast <float> (rand()/ static_cast<float> (RAND_MAX))*0.01f;
     float rand2 = static_cast <float> (rand()/ static_cast<float> (RAND_MAX))*0.01f;
     float rand3 = static_cast <float> (rand()/ static_cast<float> (RAND_MAX))*0.01f;

     CompFab::Vec3 spoint(point[0]+rand1,point[1]+rand2,point[2]+rand3);
     voxelsIntersect(voxelRes/2, voxelRes/2, voxelRes/2, spoint, !add);
     }
   }
   updateTriangleVBOs();
}
void processUpdates(){
    //cout << "Processing updates...\n";
    //cout << "Shadow pixels empty? " << shadowPixels.empty() << "\n";
    if(!shadowPixels.empty()){
        cout << "NOT EMPTY\n";    
    }
    while(!shadowPixels.empty()){
        updateLamp(shadowPixels.pop());
    }
}
void drawFilledCircle(float cx, float cy, float r){
    if (!isValidPoint(drawX, drawY)){
        return;
    }
    //fill a circle using a triangle fan
    glBegin(GL_TRIANGLE_FAN);
        glVertex2f(cx,cy);
        for (int i = 0; i <=c_res; i++)
        {
            glVertex2f(cx + 2.0*r*cos(i*c_stride*M_PI/180), cy + 2.0*r*sin(i*c_stride*M_PI/180));
        }
    glEnd();
}
/* 
  Callback that actually draws the mesh data 
*/
void displayCB()
{
    // Send off points to be processed
    processUpdates();

    // clear buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    setCamera(room_dim*0.5f, room_dim*0.7f, room_dim*2.0f, room_dim*0.5f, room_dim*0.5f,  room_dim*0.5f);
    // save the initial ModelView matrix before modifying ModelView matrix
    glPushMatrix();
    // tramsform camera
    glTranslatef(0, 0, cameraDistance);
    glRotatef(cameraAngleX, 1, 0, 0);   // pitch
    glRotatef(cameraAngleY, 0, 1, 0);   // heading

    // Enable this for mesh drawing
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    initLights();
    glColor4f(1.0f, 1.0f, 1.0f, 0.1f);
    // Draw the lamp////////////////////////////////////////////////////////////
    // Set vertex data
    glColor3f(0.6f, 0.98f, 0.95f); 
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    
    // Set normal data
    glEnableClientState(GL_NORMAL_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_normals);
    glNormalPointer(GL_FLOAT, 0, 0);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);
    glDrawElements(GL_TRIANGLES, numAct*3, GL_UNSIGNED_SHORT, 0);

    //glColor4f(0.95f, 0.0f, 0.0f, 1.0f);

    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements_test);
    //glDrawElements(GL_TRIANGLES, numInact*3, GL_UNSIGNED_SHORT, 0);

    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    // Draw the cubic room//////////////////////////////////////////////////////
    room();
  
    // Draw points
    glColor3f(1.0f, 1.0f, 1.0f); 
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    std::set<Vector3f>::iterator it;
    Vector3f point;
    for(it = pointsToDraw.begin(); it != pointsToDraw.end(); ++it){
      point = *it;
      glVertex3f(point.x(), point.y(), point.z());
    }
    glEnd();
    
    glPopMatrix();
    glutSwapBuffers();
}

void displayDrawCB()
{
    //if (!cleared){
    //    cleared = true;
    //    glClear(GL_COLOR_BUFFER_BIT);
    //}
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho (0, SCREEN_WIDTH, SCREEN_HEIGHT, 0, 0, 1);
    glMatrixMode (GL_MODELVIEW);
    glDisable(GL_DEPTH_TEST);
    
     
    glColor3f(0.7f, 0.7f, 0.7f); 
    glBegin(GL_QUADS);
    // Top left
    glVertex2f(0,0);
    glVertex2f(SCREEN_WIDTH/3,0);
    glVertex2f(SCREEN_WIDTH/3,SCREEN_HEIGHT/3);
    glVertex2f(0, SCREEN_HEIGHT/3);
    // Top right
    glVertex2f(2*SCREEN_WIDTH/3,0);
    glVertex2f(SCREEN_WIDTH,0);
    glVertex2f(SCREEN_WIDTH,SCREEN_HEIGHT/3);
    glVertex2f(2*SCREEN_WIDTH/3, SCREEN_HEIGHT/3);
    // Bottom right
    glVertex2f(2*SCREEN_WIDTH/3,2*SCREEN_HEIGHT/3);
    glVertex2f(SCREEN_WIDTH,2*SCREEN_HEIGHT/3);
    glVertex2f(SCREEN_WIDTH,SCREEN_HEIGHT);
    glVertex2f(2*SCREEN_WIDTH/3, SCREEN_HEIGHT);
    // Bottom left
    glVertex2f(0,2*SCREEN_HEIGHT/3);
    glVertex2f(SCREEN_WIDTH/3,2*SCREEN_HEIGHT/3);
    glVertex2f(SCREEN_WIDTH/3,SCREEN_HEIGHT);
    glVertex2f(0, SCREEN_HEIGHT);
    glEnd();
    // Center box
    glBegin(GL_LINE_STRIP);
    glVertex2f(SCREEN_WIDTH/3,SCREEN_HEIGHT/3);
    glVertex2f(2*SCREEN_WIDTH/3,SCREEN_HEIGHT/3);
    glVertex2f(2*SCREEN_WIDTH/3,2*SCREEN_HEIGHT/3);
    glVertex2f(SCREEN_WIDTH/3, 2*SCREEN_HEIGHT/3);
    glVertex2f(SCREEN_WIDTH/3,SCREEN_HEIGHT/3);
    glEnd();

    // User drawn shadow
    if(mouseDrawLeftDown){
      glColor3f(1.0f, 1.0f, 1.0f); 
    } else if(mouseDrawRightDown){
      glColor3f(0.0f, 0.0f, 0.0f); 
    } else {
      return;
    }
    
    drawFilledCircle(drawX, drawY, brushWidth/2);
    //glBegin(GL_POINTS);
    //for( int i = 0; i < brushWidth; i++){
    //  for( int j = 0; j < brushWidth; j++){
    //    if (isValidPoint(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j))
    //        glVertex2f(drawX-(brushWidth/2)+i, drawY-(brushWidth/2)+j);
    //  }
    //}

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
    mainWindow = initGLUT(argc, argv);
    initGL();

    // init GLUT Draw window
    drawWindow = initGLUTDraw(argc, argv);
    
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

