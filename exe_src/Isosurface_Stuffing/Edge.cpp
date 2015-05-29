
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

#include "Edge.h"
#include "Vertex.h"

//-----------------------------------------------------
Edge::Edge()
:  m_vertex1 (NULL),
   m_vertex2 (NULL),
   m_cut_point (NULL),
   m_color (BLACK_EDGE),
   m_flag (false)
{
   // Constructor
}

//-----------------------------------------------------
Edge::~Edge()
{
   // Destructor
   
   // Destroy the cut point (if it exists)
   if ( m_cut_point != NULL ) 
   {
      delete m_cut_point;
      m_cut_point = NULL;
   }
}

//-----------------------------------------------------
Vertex*
Edge::getVertex1() const
{
   // Get endpoint 1
   
   return m_vertex1;
}
   
//-----------------------------------------------------
void
Edge::setVertex1( Vertex* vertex )
{
   // Set endpoint 1
   
   m_vertex1 = vertex;
}

//-----------------------------------------------------
Vertex*
Edge::getVertex2() const
{
   // Get endpoint 2
   
   return m_vertex2;
}
   
//-----------------------------------------------------
void
Edge::setVertex2( Vertex* vertex )
{
   // Set endpoint 2
   
   m_vertex2 = vertex;
}

//-----------------------------------------------------
edge_color_type
Edge::getColor() const
{
   // Return the edge color
   
   return m_color;
}
   
//-----------------------------------------------------
void
Edge::setColor( edge_color_type color )
{
   // Set the edge color
   
   m_color = color;
}

//-----------------------------------------------------
Vertex*
Edge::getCutPoint() const
{
   // Get cut point
   
   return m_cut_point;
}
   
//-----------------------------------------------------
void
Edge::setCutPoint( Vertex* cut_point )
{
   // Set cut point
   
   m_cut_point = cut_point;
}

//-----------------------------------------------------
bool
Edge::getFlag() const
{
   return m_flag;
}
   
//-----------------------------------------------------
void
Edge::setFlag( bool flag )
{
   m_flag = flag;
}
