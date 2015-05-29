
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

#include "Octant.h"

//-----------------------------------------------------
Octant::Octant()
:  m_depth (-1),
   m_parent (NULL)
{
   // Constructor
   
   m_children.resize( 8 );
   m_children[0] = NULL;
   m_children[1] = NULL;
   m_children[2] = NULL;
   m_children[3] = NULL;
   m_children[4] = NULL;
   m_children[5] = NULL;
   m_children[6] = NULL;
   m_children[7] = NULL;
   
   m_vertices.resize( 9 );
   m_vertices[0] = NULL; 
   m_vertices[1] = NULL; 
   m_vertices[2] = NULL; 
   m_vertices[3] = NULL; 
   m_vertices[4] = NULL; 
   m_vertices[5] = NULL; 
   m_vertices[6] = NULL; 
   m_vertices[7] = NULL; 
   m_vertices[8] = NULL;
}

//-----------------------------------------------------
Octant::~Octant()
{
   // Destructor
}

//-----------------------------------------------------
int
Octant::getDepth() const
{
   return m_depth;
}

//-----------------------------------------------------
void
Octant::setDepth( int depth )
{
   m_depth = depth;
}

//-----------------------------------------------------
const VectorObj&
Octant::getGridCoords() const
{
   return m_grid_coords;
}

//-----------------------------------------------------
void
Octant::setGridCoords( const VectorObj& grid_coords )
{
   m_grid_coords = grid_coords;
}

//-----------------------------------------------------
Octant*
Octant::getParent() const
{
   return m_parent;
}

//-----------------------------------------------------
void
Octant::setParent( Octant* octant )
{
   m_parent = octant;
} 

//-----------------------------------------------------
OctantVector&
Octant::getChildren()
{
   return m_children;
}

//-----------------------------------------------------
Octant*
Octant::getChild( int num ) const
{
   if ( num < 0 || num > 7 ) {
      return NULL;
   }

   return m_children[num];
}

//-----------------------------------------------------
void
Octant::setChild( int num, Octant* octant )
{
   if ( num < 0 || num > 7 ) {
      return;
   }
   
   m_children[num] = octant;   
}

//-----------------------------------------------------
VertexVector&
Octant::getVertices()
{
   return m_vertices;
}

//-----------------------------------------------------
Vertex*
Octant::getVertex( int num ) const
{
   if ( num < 0 || num > 8 ) {
      return NULL;
   }
   
   return m_vertices[num];
}

//-----------------------------------------------------
void
Octant::setVertex( int num, Vertex* vertex )
{
   if ( num < 0 || num > 8 ) {
      return;
   }
   
   m_vertices[num] = vertex;
}

//-----------------------------------------------------
bool
Octant::containsVertex( Vertex* vertex ) const
{
   // Return if the given vertex is one of the 9 vertices of this octant.
   
   for ( int i = 0; i < 9; ++i )
   {
      if ( vertex == m_vertices[i] ) {
         return true;
      }
   }
   
   return false;
} 

//-----------------------------------------------------
bool
Octant::hasAllChildren() const
{
   // Return if this octant has all 8 children
   return
      m_children[0] != NULL &&
      m_children[1] != NULL &&
      m_children[2] != NULL &&
      m_children[3] != NULL &&
      m_children[4] != NULL &&
      m_children[5] != NULL &&
      m_children[6] != NULL &&
      m_children[7] != NULL;
}
