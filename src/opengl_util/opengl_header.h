#ifndef OPENGL_HEADER
#define OPENGL_HEADER

#if defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

#ifdef __APPLE__
#  include <GL/glew.h>
#  include <OpenGL/glu.h>
//#define USE_FREE_GLUT
#  ifdef USE_FREE_GLUT
#    include <GL/freeglut.h>
#  else
#    include <GLUT/glut.h>
#  endif
#else
#  include <GL/glew.h>
#  include <GL/glu.h>
#  include <GL/glut.h>
#endif

#endif // OPENGL_HEADER

