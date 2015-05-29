
/**
 ** Implementation of the paper:
 ** Isosurface Stuffing: Fast Tetrahedral Meshes with Good Dihedral Angles
 ** www.cs.berkeley.edu/~jrs/papers/stuffing.pdf
 ** Francois Labelle, Jonathan Richard Shewchuk
 **
 ** Implementation by Matt Gong
 ** matt.gong@gmail.com
 **
 ** ----------------------------------------------------
 **
 ** This code comes with no guarantees and must be used at the user's own
 ** risk.  I am granting permission to use this code for any purposes, but
 ** this header must remain in the code.
 **
 ** This is a visual studio 2005 project and solution.
 **
 **/

#include <unordered_map>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include <list>
#include <set>
#include <vector>
#include <string>
#include <stack>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "lib/print_macro.h"
#include "Util/Timer.h"
#include "Common.h"
#include "Tetrahedron.h"
#include "Mesh.h"
#include "Util/Util.h"
#include "Vertex.h"
#include "Edge.h"
#include "CIsoSurface.h"
#include "net_level_set.h"


// Defines

#define WINDOW_TITLE	      "Isosurface Stuffing Algorithm"
#define PI	               3.14159265
#define AXES_SIZE	         50                // Defines the world [-AXES_SIZE, +AXES_SIZE] in X and Y
//const float kScale = AXES_SIZE;
const float kScale = 1;
#define ESCAPE_KEY         0x1b              // Escape Key
#define WIN_LEFT	         30                // Lower left corner of the window
#define WIN_TOP		      30
#define WINDOW_SIZE	      700               // Initial window size
#define MIN_SCALE	         0.01              // Minimum scale factor allowed
#define ANGLE_FACTOR		   1.0               // Multiplcation factors for input interaction
#define SCALE_FACTOR		   0.005             // (These are known from previous experience)
#define LEFT_BUTTON        4                 // Active mouse buttons (or them together)
#define MIDDLE_BUTTON      2
#define RIGHT_BUTTON       1
#define BACKGROUND_COLOR	1.0, 1.0, 1.0, 0.0
#define AXES_COLOR	      0.0, 0.0, 0.0
#define AXES_WIDTH	      2.0

#define OCTREE_DEPTH 5                       // Depth of the octree
#define OCTREE_LEAF_WIDTH 0.50                // Dimension of leaf octants (approx. size of tetrahedra)
#define OCTREE_ORIGIN 0.0, 0.0, 0.0          // Origin of octree
#define LIGHT_DIRECTION 1.0, 1.0, 1.0, 0.0   // Position of directional light                                             
#define DEFAULT_SCALE 2.72                   // Default camera scale, rotation factors
#define DEFAULT_ROTX 67.0
#define DEFAULT_ROTZ 154.0
#define CAMERA_TARGET 0.0, 0.0, 0.0          // Camera target

#define SPHERE_RADIUS 20                     // Test surface sphere radius 
#define SURFACE_OCTANT 20.0, 0.0, 0.0        // 3D coordinate for an octant on surface (initial guess)


// Enumerations

enum transformation_mode_type {                 // Left button - rotating or scaling
  ROTATE_MODE = 0,
  SCALE_MODE
};


void WriteTetGen(const char *file_name, std::vector<double> &verts, std::vector<int> &tets) {
  char filename[1024];
  sprintf(filename, "%s.node", file_name);
  FILE *fp = fopen(filename, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", filename);
    return;
  }
  ASSERT(verts.size() % 3 == 0);
  ASSERT(tets.size() % 4 == 0);
  int vertex_num = int(verts.size()) / 3;
  int tet_number = int(tets.size()) / 4;
  fprintf(fp, "%d %d %d %d\n", vertex_num, 3, 0, 0);
  for (int i = 0; i < vertex_num; i++)
    fprintf(fp, "%d %f %f %f\n", i + 1, verts[i * 3 + 0], verts[i * 3 + 1], verts[i * 3 + 2]);
  fclose(fp);

  sprintf(filename, "%s.ele", file_name);
  fp = fopen(filename, "w+");
  if (fp == NULL)	{
    printf("ERROR: file %s not open.\n", filename);
    return;
  }
  fprintf(fp, "%d %d %d\n", tet_number, 4, 0);

  for (int i = 0; i < tet_number; i++)
    fprintf(fp, "%d %d %d %d %d\n", i + 1, tets[i * 4 + 0] + 1, tets[i * 4 + 1] + 1, tets[i * 4 + 2] + 1, tets[i * 4 + 3] + 1);
  fclose(fp);
}
// Global Variables

int	      window_id = 0;	                     // Window id for the top level window

int	      axes_list = 0;		                  // Axes display list
bool	      axes_visible = false;               // Axes on/off

int         octree_list = 0;                    // Octree (leaves only) display list
bool        octree_visible = false;

int         octree_parents_list = 0;            // Octree parent octants display list
bool        octree_parents_visible = true;

int         cut_point_list = 0;                 // The bcc lattice cut point display list
bool        cut_point_visible = false;

int         mesh_list = 0;                      // The mesh display list
bool        mesh_visible = true;

int         wire_frame_list = 0;                // The wireframe mesh display list
bool        wire_frame_visible = true;

transformation_mode_type transformation_mode = ROTATE_MODE;    // Rotate or Scale mode

int         active_button = 0;		            // Current button that is down
double	   scale_factor = 0.0;                 // Scene scale factor
double      rotation_x = 0.0;                   // Rotation angle in degrees in X
double      rotation_z = 0.0;                   // Rotation angle in degrees in Z
int	      mouse_x = 0;                        // Mouse coordinates in X
int         mouse_y = 0;                        // Mouse coordinates in Y

Mesh        mesh;                               // Tetrahedral mesh


// Function prototypes

void animate();
void display();
void initGraphics();
void initDisplayLists();
void keyboard( unsigned char, int, int );
void mouseButton( int button, int state, int x, int y );
void mouseMotion( int x, int y );
void quit();
void reset();
void resize( int width, int height );
void visibility( int state );
void drawAxes( double length );
void rasterStringBig( double x, double y, double z, char *s );
void rasterStringSmall( double x, double y, double z, char *s );
void strokeString( double x, double y, double z, double ht, char *s );

//-------------------------------------------------------------------------
double
phi( const VectorObj& coords ) {
  // The cut function, or signed distance function, that specifies
  // the isosurface.
  // Postive values are inside the surface, negative values are outside,
  // and 0 is on the surface.

  // Test sphere
  return SPHERE_RADIUS - sqrt(coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2]) + 2;
}

//-------------------------------------------------------------------------
int
main( int argc, char *argv[] ) {

  // Turn on the glut package:
  // (do this before checking argc and argv since it might
  // pull some command line arguments out)
  glutInit( &argc, argv );

  // Read the command line
  for ( int i = 1; i < argc; ++i ) {
    fprintf( stderr, "Unknown argument: '%s'\n", argv[i] );
    fprintf( stderr, "Usage: %s [-D]\n", argv[0] );
  }

  // Print keys
  std::cout << "A = Axes\n";
  std::cout << "Q = Quit\n";
  std::cout << "R = Rotate\n";
  std::cout << "S = Scale\n";
  std::cout << "E = Reset\n";
  std::cout << "L = Octree\n";
  std::cout << "V = Octree Leaves\n";
  std::cout << "C = Cut Points\n";
  std::cout << "M = Mesh\n";
  std::cout << "W = Wireframe Mesh\n";

  // Setup all graphics - including callbacks
  initGraphics();

  //  const float kSize = AXES_SIZE / 8.0;
  const float kSize = 1.0;
  double width = 1.0;
  double height = 2.4;
  double radius = 0.03;
  int nx = 5;
  int ny = 12;
  double x0 = -width / 2 + 0.11;
  double y0 = -height / 2 + 0.1;
  double dx = 2 * radius + 0.135;
  double dy = 2 * radius + 0.135;
  NetLevelSet net(width, height, radius, nx, ny, x0, y0, dx, dy);
  //  NetLevelSet net(30, 30, 10, 1, 1, 0, 0, 0, 0);
  // Perform the meshing algorithm
  mesh.generateMesh(
    //      phi,
    net,
    5, 0.061,
    VectorObj(OCTREE_ORIGIN),
    VectorObj(1.0 * kSize, 0, 0));

  // Get the tetrahedra and print size
  const TetrahedronVector& tetrahedra = mesh.getTetrahedra();
  std::cout << "Mesh contains " << (long) tetrahedra.size() << " tetrahedra." << std::endl;
    if (0) {
    TetrahedronVector& tet = mesh.m_tetrahedra;
    VertexVector& v = mesh.m_vertices;
    std::vector<double> verts;
    std::vector<int> tets;
    tets.resize(tet.size() * 4);
    verts.resize(v.size() * 3);
    std::unordered_map<Vertex*, int> vert_idx_map;
    for (int i = 0; i < int(v.size()); ++i) {
      verts[i * 3 + 0] = v[i]->x();
      verts[i * 3 + 1] = v[i]->y();
      verts[i * 3 + 2] = v[i]->z();
      vert_idx_map[v[i]] = i;
    }
    for (int t = 0; t < int(tet.size()); ++t) {
      Vertex* pos[4];
      tet[t]->getVertices(pos[0], pos[1], pos[2], pos[3]);
      tets[t * 4 + 0] = vert_idx_map[pos[0]];
      tets[t * 4 + 1] = vert_idx_map[pos[1]];
      tets[t * 4 + 2] = vert_idx_map[pos[2]];
      tets[t * 4 + 3] = vert_idx_map[pos[3]];
    }
    WriteTetGen("/Users/dj/mesh", verts, tets);
    exit(0);
  }

  // Create the display lists
  initDisplayLists();

  // Initialize the transformations and colors
  // (and post a redisplay)
  reset();

  // draw the scene and wait for interaction
  // (will never return)
  glutMainLoop();

  return 0;
}

//-------------------------------------------------------------------------
void
display() {
#if 0
#else
  // Display the scene

  // Set the current window
  glutSetWindow( window_id );

  // Clear background
  glDrawBuffer( GL_BACK );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glEnable( GL_DEPTH_TEST );

  // Set ambient lighting
  GLfloat ambient_global[] = { 0.4, 0.4, 0.4, 1.0 };
  glLightModelfv( GL_LIGHT_MODEL_AMBIENT, ambient_global );

  // Create a light
  glEnable( GL_LIGHT0 );
  GLfloat diffuse0[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat specular0[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat position0[] = { LIGHT_DIRECTION };
  glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse0 );
  glLightfv( GL_LIGHT0, GL_SPECULAR, specular0 );
  glLightfv( GL_LIGHT0, GL_POSITION, position0 );

  // Set the shading model to smooth
  glShadeModel( GL_SMOOTH );

  // Normalize all surface normals so our mouse zooming doesn't change the lighting
  glEnable( GL_NORMALIZE );

  // Set the viewport to a square centered in the window
  int window_width = glutGet( GLUT_WINDOW_WIDTH );
  int window_height = glutGet( GLUT_WINDOW_HEIGHT );
  int minimum_dimension = window_width < window_height ? window_width : window_height;
  int viewport_ll_x = ( window_width - minimum_dimension ) / 2;
  int viewport_ll_y = ( window_height - minimum_dimension ) / 2;
  glViewport( viewport_ll_x, viewport_ll_y, minimum_dimension, minimum_dimension );

  // Setup the projection matrix:
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( 80.0, 1.0,	0.1, 400.0 );

  // Setup modelview matrix:
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  // Camera is offset looking at the octree origin
  gluLookAt( 50.0, 50.0, 50.0, CAMERA_TARGET, 0.0, 0.0, 1.0 );

  // Rotate and scale scene according to user interaction
  VectorObj camera_target( CAMERA_TARGET );
  glTranslatef( camera_target[0], camera_target[1], camera_target[2] );
  glRotatef( rotation_z, 0.0, 0.0, 1.0 );
  glRotatef( rotation_x, 1.0, 0.0, 0.0 );
  glScalef( scale_factor, scale_factor, scale_factor );
  glTranslatef( -camera_target[0], -camera_target[1], -camera_target[2] );

  // Enable lighting
  glEnable( GL_LIGHTING );

  // Draw mesh
  if ( mesh_visible ) {
    glCallList( mesh_list );
  }

  // Disable lighting
  glDisable( GL_LIGHTING );

  // Enable anti-aliased lines
  glEnable( GL_LINE_SMOOTH );
  glEnable( GL_BLEND );
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  glHint( GL_LINE_SMOOTH_HINT, GL_DONT_CARE );

  // Draw the axes
  if ( axes_visible ) {
    glCallList( axes_list );
  }
  glPushMatrix();
  glScalef(kScale, kScale, kScale);
  glEnable(GL_NORMALIZE);

  // Draw the Octree
  if ( octree_visible ) {
    glCallList( octree_list );

    if ( octree_parents_visible ) {
      glCallList( octree_parents_list );
    }
  }

  // Draw the cut points
  if ( cut_point_visible ) {
    glCallList( cut_point_list );
  }

  // Draw the wireframe mesh
  if ( wire_frame_visible ) {
    // Enable a polygon offset to remove ugly "stiching" artificats
    // when drawing the wireframe over the mesh.
    glEnable( GL_POLYGON_OFFSET_LINE );
    glPolygonOffset( 0.0, -20.0 );

    glCallList( wire_frame_list );

    glDisable(GL_POLYGON_OFFSET_FILL);
  }
  glPopMatrix();
  // Swap the double-buffered framebuffers
  glutSwapBuffers();

  // Flush the scene
  glFlush();
#endif
}

//-------------------------------------------------------------------------
void
initGraphics() {
  // Initialize the glut and OpenGL libraries.
  // Setup display lists and callback functions.

  // Setup the display mode
  // ( must be done before glutCreateWindow() )
  glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );

  // Set the initial window configuration
  glutInitWindowSize( WINDOW_SIZE, WINDOW_SIZE );
  glutInitWindowPosition( WIN_LEFT, WIN_TOP );

  // Open the window and set its title
  window_id = glutCreateWindow( WINDOW_TITLE );
  glutSetWindowTitle( WINDOW_TITLE );

  // Setup the clear color
  glClearColor( BACKGROUND_COLOR );

  /* Setup callback functions:
    *
   * DisplayFunc -- redraw the window
   * ReshapeFunc -- handle the user resizing the window
   * KeyboardFunc -- handle a keyboard input
   * MouseFunc -- handle the mouse button going down or up
   * MotionFunc -- handle the mouse moving with a button down
   * PassiveMotionFunc -- handle the mouse moving with a button up
   * VisibilityFunc -- handle a change in window visibility
   * EntryFunc	-- handle the cursor entering or leaving the window
   * SpecialFunc -- handle special keys on the keyboard
   * SpaceballMotionFunc -- handle spaceball translation
   * SpaceballRotateFunc -- handle spaceball rotation
   * SpaceballButtonFunc -- handle spaceball button hits
   * ButtonBoxFunc -- handle button box hits
   * DialsFunc -- handle dial rotations
   * TabletMotionFunc -- handle digitizing tablet motion
   * TabletButtonFunc -- handle digitizing tablet button hits
   * MenuStateFunc -- declare when a pop-up menu is in use
   * IdleFunc -- what to do when nothing else is going on
   * TimerFunc -- trigger something to happen every so often
    *
    */
  glutSetWindow( window_id );
  glutDisplayFunc( display );
  glutReshapeFunc( resize );
  glutKeyboardFunc( keyboard );
  glutMouseFunc( mouseButton );
  glutMotionFunc( mouseMotion );
  glutVisibilityFunc( visibility );
  //  glutPassiveMotionFunc( NULL );
  //  glutEntryFunc( NULL );
  //	glutSpecialFunc( NULL );
  //	glutSpaceballMotionFunc( NULL );
  //	glutSpaceballRotateFunc( NULL );
  //	glutSpaceballButtonFunc( NULL );
  //	glutButtonBoxFunc( NULL );
  //	glutDialsFunc( NULL );
  //	glutTabletMotionFunc( NULL );
  //	glutTabletButtonFunc( NULL );
  //	glutMenuStateFunc( NULL );
  //	glutTimerFunc( 0, NULL, 0 );
}

//-------------------------------------------------------------------------
void
initDisplayLists() {
  // Initialize the display lists

  // Axes
  axes_list = glGenLists( 1 );
  glNewList( axes_list, GL_COMPILE );
  glColor3f( AXES_COLOR );
  glLineWidth( AXES_WIDTH );
  drawAxes( AXES_SIZE );
  glEndList();

  // Octree
  octree_list = glGenLists( 1 );
  glNewList( octree_list, GL_COMPILE );
  mesh.drawOctreeLeaves();
  glEndList();

  // Octree parents
  octree_parents_list = glGenLists( 1 );
  glNewList( octree_parents_list, GL_COMPILE );
  mesh.drawOctreeParents();
  glEndList();

  // Cut Points
  cut_point_list = glGenLists( 1 );
  glNewList( cut_point_list, GL_COMPILE );
  mesh.drawCutPoints();
  glEndList();

  // Mesh
  mesh_list = glGenLists( 1 );
  glNewList( mesh_list, GL_COMPILE );
  mesh.drawMesh();
  glEndList();

  // Wireframe Mesh
  wire_frame_list = glGenLists( 1 );
  glNewList( wire_frame_list, GL_COMPILE );
  mesh.drawWireframeMesh();
  glEndList();
}

//-------------------------------------------------------------------------
void
keyboard( unsigned char c, int x, int y ) {
  // Keyboard keys

  switch ( c ) {
    case 'a':
    case 'A':
      axes_visible = ! axes_visible;
      break;

    case 'q':
    case 'Q':
    case ESCAPE_KEY:
      quit();
      break;

    case 'r':
    case 'R':
      transformation_mode = ROTATE_MODE;
      break;

    case 's':
    case 'S':
      transformation_mode = SCALE_MODE;
      break;

    case 'e':
    case 'E':
      reset();
      break;

    case 'l':
    case 'L':
      octree_visible = ! octree_visible;
      break;

    case 'v':
    case 'V':
      octree_parents_visible = ! octree_parents_visible;
      break;

    case 'c':
    case 'C':
      cut_point_visible = ! cut_point_visible;
      break;

    case 'm':
    case 'M':
      mesh_visible = ! mesh_visible;
      break;

    case 'w':
    case 'W':
      wire_frame_visible = ! wire_frame_visible;
      break;

    default:
      break;
  }

  glutSetWindow( window_id );
  glutPostRedisplay();
}

//-------------------------------------------------------------------------
void
mouseButton( int button, int state, int x, int y ) {
  // User pressed mouse button

  // Get the proper button bit mask
  int bit_mask = 0;
  switch ( button ) {
    case GLUT_LEFT_BUTTON:
      bit_mask = LEFT_BUTTON;
      break;

    case GLUT_MIDDLE_BUTTON:
      bit_mask = MIDDLE_BUTTON;
      break;

    case GLUT_RIGHT_BUTTON:
      bit_mask = RIGHT_BUTTON;
      break;

    default:
      bit_mask = 0;
      fprintf( stderr, "Unknown mouse button: %d\n", button );
      break;
  }

  // Set or clear the proper bit
  if ( state == GLUT_DOWN ) {
    mouse_x = x;
    mouse_y = y;
    active_button |= bit_mask;
  } else {
    active_button &= ~bit_mask;
  }
}

//-------------------------------------------------------------------------
void
mouseMotion( int x, int y ) {
  // User moved mouse

  // Delta mouse coordinates
  int delta_x = x - mouse_x;
  int delta_y = y - mouse_y;

  // Left mouse button
  if ( active_button & LEFT_BUTTON ) {
    if ( transformation_mode == ROTATE_MODE ) {
      rotation_x += ( ANGLE_FACTOR * delta_y );
      rotation_z += ( ANGLE_FACTOR * delta_x );
    } else {
      scale_factor += SCALE_FACTOR * (double)( delta_x - delta_y );

      // Do not invert
      if ( scale_factor < MIN_SCALE ) {
        scale_factor = MIN_SCALE;
      }
    }
  }

  // Middle mouse button
  if ( active_button & MIDDLE_BUTTON ) {
    scale_factor += SCALE_FACTOR * (double)( delta_x - delta_y );

    // Do not invert
    if ( scale_factor < MIN_SCALE ) {
      scale_factor = MIN_SCALE;
    }
  }

  // Set new current position
  mouse_x = x;
  mouse_y = y;

  glutSetWindow( window_id );
  glutPostRedisplay();
}

//-------------------------------------------------------------------------
void
quit() {
  // Quit the program gracefully

  // Destroy the display lists
  glDeleteLists( axes_list, 1 );
  glDeleteLists( octree_list, 1 );
  glDeleteLists( octree_parents_list, 1 );
  glDeleteLists( mesh_list, 1 );
  glDeleteLists( wire_frame_list, 1 );
  glDeleteLists( cut_point_list, 1 );

  glFinish();

  // Destroy the window
  glutDestroyWindow( window_id );

  exit( 0 );
}

//-------------------------------------------------------------------------
void
reset() {
  // Reset the scene

  active_button = 0;

  // Set default rotation/scale
  scale_factor = DEFAULT_SCALE;
  rotation_x = DEFAULT_ROTX;
  rotation_z = DEFAULT_ROTZ;

  glutSetWindow( window_id );
  glutPostRedisplay();
}

//-------------------------------------------------------------------------
void
resize( int width, int height ) {
  // Window resized

  glutSetWindow( window_id );
  glutPostRedisplay();
}

//-------------------------------------------------------------------------
void
visibility( int state ) {
  // Window visibility change

  if ( state == GLUT_VISIBLE ) {
    glutSetWindow( window_id );
    glutPostRedisplay();
  }
}

//-------------------------------------------------------------------------
void
drawAxes( double length ) {

  // The stroke characters 'X' 'Y' 'Z'
  double xx[] = { 0.0, 1.0, 0.0, 1.0 };
  double xy[] = { -0.5, 0.5, 0.5, -0.5 };
  int xorder[] = { 1, 2, -3, 4 };
  double yx[] = { 0.0, 0.0, -0.5, 0.5 };
  double yy[] = { 0.0, 0.6, 1.0, 1.0 };
  int yorder[] = { 1, 2, 3, -2, 4 };
  double zx[] = { 1.0, 0.0, 1.0, 0.0, 0.25, 0.75 };
  double zy[] = { 0.5, 0.5, -0.5, -0.5, 0.0, 0.0 };
  int zorder[] = { 1, 2, 3, 4, -5, 6 };

  // Fraction of the length to use as height of the characters
  double length_fraction = 0.10;

  // Fraction of length to use as start location of the characters
  double base_fraction = 1.10;

  glBegin( GL_LINE_STRIP );
  glVertex3f( length, 0., 0. );
  glVertex3f( 0., 0., 0. );
  glVertex3f( 0., length, 0. );
  glEnd();
  glBegin( GL_LINE_STRIP );
  glVertex3f( 0., 0., 0. );
  glVertex3f( 0., 0., length );
  glEnd();

  // Character scale factor and start location
  double fact = length_fraction * length;
  double base = base_fraction * length;

  int j = 0;
  glBegin( GL_LINE_STRIP );
  for ( int i = 0; i < 4; i++ ) {
    j = xorder[i];
    if ( j < 0 ) {

      glEnd();
      glBegin( GL_LINE_STRIP );
      j = -j;
    }
    j--;
    glVertex3f( base + fact * xx[j], fact * xy[j], 0.0 );
  }
  glEnd();

  glBegin( GL_LINE_STRIP );
  for ( int i = 0; i < 5; i++ ) {
    j = yorder[i];
    if ( j < 0 ) {

      glEnd();
      glBegin( GL_LINE_STRIP );
      j = -j;
    }
    j--;
    glVertex3f( fact * yx[j], base + fact * yy[j], 0.0 );
  }
  glEnd();

  glBegin( GL_LINE_STRIP );
  for ( int i = 0; i < 6; i++ ) {
    j = zorder[i];
    if ( j < 0 ) {

      glEnd();
      glBegin( GL_LINE_STRIP );
      j = -j;
    }
    j--;
    glVertex3f( 0.0, fact * zy[j], base + fact * zx[j] );
  }
  glEnd();
}

//-------------------------------------------------------------------------
void
rasterStringBig( double x, double y, double z, double *s ) {
  // Use glut to display a string of characters using a raster font

  // One character to print
  char c;

  glRasterPos3f( x, y, z );
  for ( ; ( c = *s ) != '\0'; s++ ) {
    glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, c );
  }
}

//-------------------------------------------------------------------------
void
rasterStringSmall( double x, double y, double z, double *s ) {
  char c;			/* one character to print		*/

  glRasterPos3f( x, y, z );
  for ( ; ( c = *s ) != '\0'; s++ ) {
    glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, c );
  }
}

//-------------------------------------------------------------------------
void
strokeString( double x, double y, double z, double ht, char *s ) {
  // Use glut to display a string of characters using a stroke font

  // Character to print
  char c;

  // Scale factor
  double sf;

  glPushMatrix();
  glTranslatef( x, y, z );
  sf = ht / ( 119.05 + 33.33 );
  glScalef( sf, sf, sf );
  for ( ; ( c = *s ) != '\0'; s++ ) {
    glutStrokeCharacter( GLUT_STROKE_ROMAN, c );
  }
  glPopMatrix();
}
