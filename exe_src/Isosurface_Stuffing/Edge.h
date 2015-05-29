
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

#ifndef _EDGE_H_
#define _EDGE_H_

#include <list>
#include "Common.h"

// Forward Class Declarations
class Vertex;

// Edge Class
class Edge
{
public:
  
   Edge();
   virtual ~Edge();
   
   Vertex*
   getVertex1() const; 
   void
   setVertex1( Vertex* vertex );

   Vertex*
   getVertex2() const; 
   void
   setVertex2( Vertex* vertex );

   edge_color_type
   getColor() const; 
   void
   setColor( edge_color_type color );

   Vertex*
   getCutPoint() const; 
   void
   setCutPoint( Vertex* cut_point );
   
   bool
   getFlag() const;
   void
   setFlag( bool flag );

private:
   
   Vertex* m_vertex1;               // Two end vertices
   Vertex* m_vertex2;
   edge_color_type m_color;         // Edge Color: Black = axis aligned, Red = diagonal, Blue = split edge   
   Vertex* m_cut_point;             // Optional "cut" point - NULL if it doesn't exist
                                    // This memory is managed by this edge.
   bool m_flag;                     // Mark this edge as it is processed
      
};

#endif 
