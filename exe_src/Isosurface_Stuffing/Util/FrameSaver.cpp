#ifdef WIN32
#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif


#include "FrameSaver.h"

FrameSaver::FrameSaver()
{
	m_record = 0 ;
	m_pixels = NULL ;
	m_frameCount = 0 ;
}


void FrameSaver::Toggle(int w)
{
	if( m_record == 0 ) 
		StartRecord(w) ;
	else 
		m_record = 0 ;
}

void FrameSaver::StartRecord(int width)
{
	m_record = 1 ;
	m_frameCount = 0 ;
}


int FrameSaver::DumpPPM(int width, int height)
{
	if( m_pixels ) delete [] m_pixels ;
	m_pixels = new unsigned char [3*width];
    if(	m_pixels == NULL )
    {
		fprintf(stderr,"Cannot allocate	enough memory\n") ;
		return  -1;
    }

  char fname[1024] ;
	if( m_record == 0 ) // one time save
    strcpy(fname, "scene.ppm") ;
	else				// animation recording
	{
    sprintf(fname, "frame%d.ppm", m_frameCount) ;
		m_frameCount++ ;
	}
  FILE *fp = fopen(fname, "wb") ;
  if( fp != 0 )
	{
		fprintf(stderr, "Cannot open file %s\n", fname) ;
		return -1 ;
	}
	DumpPPM(fp,width,height) ;
	fclose(fp) ;
	return 1 ;
}

// dumps a PPM raw (P6) file on an already allocated memory array
void FrameSaver::DumpPPM(FILE *fp, int width, int height)
{
    const int maxVal=255;
    register int y;
    

     //printf("width = %d height = %d\n",width, height) ;

    fprintf(fp,	"P6 ");
    fprintf(fp,	"%d %d ", width, height);
    fprintf(fp,	"%d\n",	maxVal);

   	glReadBuffer(GL_FRONT) ; 
	for	( y = height-1;	y>=0; y-- ) 
	{
		glReadPixels(0,y,width,1,GL_RGB,GL_UNSIGNED_BYTE, 
			(GLvoid *) m_pixels);
		fwrite(m_pixels, 3, width, fp);
	}
}
