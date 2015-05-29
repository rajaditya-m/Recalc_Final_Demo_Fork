
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

#include "Vertex.h"
#include "Edge.h"
#include "Tetrahedron.h"

//-----------------------------------------------------
Vertex::Vertex()
:  m_coords (0.0, 0.0, 0.0),
   m_grid_coords (0.0, 0.0, 0.0),
   m_color (BLACK_VERTEX),
   m_f (-1.0),
   m_f_computed (false),
   m_in_subset (false)
{
   // Constructor
}

//-----------------------------------------------------
Vertex::~Vertex()
{
   // Destructor
}

//-----------------------------------------------------
VectorObj&
Vertex::getCoords()
{
   return m_coords;
}

//-----------------------------------------------------
void
Vertex::setCoords( const VectorObj& coords )
{
   m_coords = coords;
}

//-----------------------------------------------------
VectorObj&
Vertex::getGridCoords()
{
   return m_grid_coords;
}

//-----------------------------------------------------
void
Vertex::setGridCoords( const VectorObj& grid_coords )
{
   m_grid_coords = grid_coords;
}

//-----------------------------------------------------
EdgeList&
Vertex::getEdgeList()
{
   // Return edge list
   
   return m_edge_list;
}
   
//-----------------------------------------------------
void
Vertex::setEdgeList( const EdgeList& edge_list )
{
   // Set the edge list
   
   m_edge_list = edge_list;
}

//-----------------------------------------------------
TetrahedronVector&
Vertex::getTetrahedronVector()
{
   return m_tetrahedron_vector;
}

//-----------------------------------------------------
void
Vertex::setTetrahedronVector( const TetrahedronVector& tetrahedron_vector )
{
   m_tetrahedron_vector = tetrahedron_vector;
}

//-----------------------------------------------------
Edge*
Vertex::getEdge( Vertex* end_vertex ) const
{
   // Get the edge that connects this vertex to the end vertex
   
   // Loop through all adjacent edges to this vertex
   Edge* edge = NULL;
   EdgeListConstIterator edge_iter = m_edge_list.begin();
   for ( ; edge_iter != m_edge_list.end(); ++edge_iter )
   {
      edge = *edge_iter;
      
      // Return if we found our edge
      if ( edge->getVertex1() == this && edge->getVertex2() == end_vertex ||
           edge->getVertex1() == end_vertex && edge->getVertex2() == this )
      {
         return edge;
      }
   }
   
   // Not found
   return NULL;   
}

//-----------------------------------------------------
vertex_color_type
Vertex::getColor() const
{
   // Return the vertex color
   
   return m_color;
}
   
//-----------------------------------------------------
void
Vertex::setColor( vertex_color_type color )
{
   // Set the vertex color
   
   m_color = color;
}

//-----------------------------------------------------
double
Vertex::getF() const
{
   // Return value of vertex f
      
   return m_f;
}
   
//-----------------------------------------------------
void 
Vertex::setF( double f )
{
   // Set value of vertex f

   m_f = f;
   
   m_f_computed = true;
}

//-----------------------------------------------------
bool
Vertex::getFComputed() const
{
   return m_f_computed;
}

//-----------------------------------------------------
bool
Vertex::getInSubset() const
{
   return m_in_subset;
}

//-----------------------------------------------------
void
Vertex::setInSubset( bool in_subset )
{
   m_in_subset = in_subset;
}