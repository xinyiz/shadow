#include <GL/glut.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>
#include <vecmath.h>
using namespace std;


// Globals


// This is the list of points (3D vectors)
vector<Vector3f> vecv;


// This is the list of normals (also 3D vectors)
vector<Vector3f> vecn;


// This is the list of faces (indices into vecv and vecn)
vector<vector<unsigned> > vecf;




// You will need more global variables to implement color and position changes
int color = 0;
float lightposx = 1.0f;
float lightposy = 1.0f;
int MAX_BUFFER_SIZE = 1000;
// These are convenience functions which allow us to call OpenGL
// methods on Vec3d objects
inline void glVertex(const Vector3f &a)
{ glVertex3fv(a); }


inline void glNormal(const Vector3f &a)
{ glNormal3fv(a); }




// This function is called whenever a "Normal" key press is received.
void keyboardFunc( unsigned char key, int x, int y )
{
	switch ( key )
	{
		case 27: // Escape key
			exit(0);
			break;
		case 'c':
			color ++;
			if (color == 4){
				color = 0;
			}
			break;
		default:
			cout << "Unhandled key press " << key << "." << endl;
	}


	// this will refresh the screen so that the user sees the color change
	glutPostRedisplay();
}


// This function is called whenever a "Special" key press is received.
// Right now, it's handling the arrow keys.
void specialFunc( int key, int x, int y )
{
	switch ( key )
	{
		case GLUT_KEY_UP:
			lightposy -= 0.5f;
			break;
		case GLUT_KEY_DOWN:
			lightposy += 0.5;
			break;
		case GLUT_KEY_LEFT:
			lightposx -= 0.5;
			break;
		case GLUT_KEY_RIGHT:
			lightposx += 0.5;
			break;
	}




	// this will refresh the screen so that the user sees the light position
	glutPostRedisplay();
}


// This function is responsible for displaying the object.
void drawScene(void)
{
	//int i;
	int a;
	int c;
	int d;
	int f;
	int g;
	int i;
	// Clear the rendering window
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	// Rotate the image
	glMatrixMode( GL_MODELVIEW );  // Current matrix affects objects positions
	glLoadIdentity();                  // Initialize to the identity


	// Position the camera at [0,0,5], looking at [0,0,0],
	// with [0,1,0] as the up direction.
	gluLookAt(0.0, 0.0, 5.0,
			0.0, 0.0, 0.0,
			0.0, 1.0, 0.0);


	// Set material properties of object


	// Here are some colors you might use - feel free to add more
	GLfloat diffColors[4][4] = { {0.5, 0.5, 0.9, 1.0},
		{0.9, 0.5, 0.5, 1.0},
		{0.5, 0.9, 0.3, 1.0},
		{0.3, 0.8, 0.9, 1.0} };

	// Here we use the first color entry as the diffuse color
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, diffColors[color]);


	// Define specular color and shininess
	GLfloat specColor[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat shininess[] = {100.0};


	// Note that the specular color and shininess can stay constant
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specColor);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);

	// Set light properties


	// Light color (RGBA)
	GLfloat Lt0diff[] = {1.0,1.0,1.0,1.0};
	// Light position
	GLfloat Lt0pos[] = {lightposx, lightposy, 5.0f, 1.0f};


	glLightfv(GL_LIGHT0, GL_DIFFUSE, Lt0diff);
	glLightfv(GL_LIGHT0, GL_POSITION, Lt0pos);
	for(unsigned int p=0; p < vecf.size(); p++) {
		vector<unsigned> &vf = vecf[p];
		a = vf[0];
		c = vf[1];
		d = vf[2];
		f = vf[3];
		g = vf[4];
		i = vf[5];
		glBegin(GL_TRIANGLES);
		glNormal3d(vecn[c-1][0],vecn[c-1][1],vecn[c-1][2]);
		glVertex3d(vecv[a-1][0],vecv[a-1][1],vecv[a-1][2]);
		glNormal3d(vecn[f-1][0],vecn[f-1][1],vecn[f-1][2]);
		glVertex3d(vecv[d-1][0],vecv[d-1][1],vecv[d-1][2]);
		glNormal3d(vecn[i-1][0],vecn[i-1][1],vecn[i-1][2]);
		glVertex3d(vecv[g-1][0],vecv[g-1][1],vecv[g-1][2]);
		glEnd();
	}
	//glutSolidTeapot(1.0);
	// Dump //the image to the screen.
	glutSwapBuffers();




}


// Initialize OpenGL's rendering modes
void initRendering()
{
	glEnable(GL_DEPTH_TEST);   // Depth testing must be turned on
	glEnable(GL_LIGHTING);         // Enable lighting calculations
	glEnable(GL_LIGHT0);           // Turn on light #0.
}


// Called when the window is resized
// w, h - width and height of the window in pixels.
void reshapeFunc(int w, int h)
{
	// Always use the largest square viewport possible
	if (w > h) {
		glViewport((w - h) / 2, 0, h, h);
	} else {
		glViewport(0, (h - w) / 2, w, w);
	}


	// Set up a perspective view, with square aspect ratio
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// 50 degree fov, uniform aspect ratio, near = 1, far = 100
	gluPerspective(50.0, 1.0, 1.0, 100.0);
}


void loadInput()
{
	Vector3f v;
	Vector3f n;
	string s;
	int index = 0;
	char buffer[MAX_BUFFER_SIZE];
	while(cin.getline(buffer, MAX_BUFFER_SIZE) != 0){
		stringstream ss(buffer);
		ss >> s;
		if (s == "v"){
			ss >> v[0] >> v[1] >> v[2];
			vecv.push_back(Vector3f(v[0],v[1],v[2]));
		}
		else if(s == "vn"){
			ss >> n[0] >> n[1] >> n[2];
			vecn.push_back(Vector3f(n[0],n[1],n[2]));
		}
		else if(s == "f"){
			string str = buffer;
			string word;
			stringstream stream(str.substr(1));
			vector<unsigned> line;
			index = 0;
			while( getline(stream, word, '/') ){
				if (index == 0 | index == 6){
					line.push_back(atoi(word.c_str()));
				}
				else if(index == 2 | index == 4){
					istringstream iss(word);
					string first;
					string last;
					iss >> first >> last;
					line.push_back(atoi(first.c_str()));
					line.push_back(atoi(last.c_str()));
				}
				index++;
			}
			vecf.push_back(line);
		}
	}


}


// Main routine.
// Set up OpenGL, define the callbacks and start the main loop
int main( int argc, char** argv )
{
	loadInput();


	glutInit(&argc,argv);


	// We're going to animate it, so double buffer
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );


	// Initial parameters for window position and size
	glutInitWindowPosition( 60, 60 );
	glutInitWindowSize( 360, 360 );
	glutCreateWindow("Assignment 0");


	// Initialize OpenGL parameters.
	initRendering();


	// Set up callback functions for key presses
	glutKeyboardFunc(keyboardFunc); // Handles "normal" ascii symbols
	glutSpecialFunc(specialFunc);   // Handles "special" keyboard keys


	// Set up the callback function for resizing windows
	glutReshapeFunc( reshapeFunc );


	// Call this whenever window needs redrawing
	glutDisplayFunc( drawScene );


	// Start the main loop.  glutMainLoop never returns.
	glutMainLoop( );


	return 0;    // This line is never reached.
}
