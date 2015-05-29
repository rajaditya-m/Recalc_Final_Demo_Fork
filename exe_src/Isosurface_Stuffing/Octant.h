
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

#ifndef _OCTANT_H_
#define _OCTANT_H_

#include "Common.h"
#include "Util/vectorObj.h"

// Forward Class Declarations
class Vertex;
class Edge;

// Octant class
class Octant
{
public:
  
   Octant();
   virtual ~Octant();

   int
   getDepth() const;
   void
   setDepth( int depth );

   const VectorObj&
   getGridCoords() const;   
   void
   setGridCoords( const VectorObj& grid_coords );

   Octant*
   getParent() const;  
   void
   setParent( Octant* octant ); 

   OctantVector&
   getChildren();

   Octant*
   getChild( int num ) const;   
   void
   setChild( int num, Octant* octant );

   VertexVector&
   getVertices();

   Vertex*
   getVertex( int num ) const; 
   void
   setVertex( int num, Vertex* vertex );
   
   bool
   containsVertex( Vertex* vertex ) const;
   
   bool
   hasAllChildren() const; 

private:
   
   int m_depth;               // Octant Depth
   VectorObj m_grid_coords;   // Grid coordinates at current depth (i.e. (0,0,0), (0,0,1), etc.) 

   Octant* m_parent;          // Parent octant

   OctantVector m_children;   // Children [0-7]  
                              // Child 0 = octant (0, 0, 0)
                              // Child 1 = octant (1, 0, 0)
                              // Child 2 = octant (0, 0, 1)
                              // Child 3 = octant (1, 0, 1)
                              // Child 4 = octant (0, 1, 0)
                              // Child 5 = octant (1, 1, 0)
                              // Child 6 = octant (0, 1, 1)
                              // Child 7 = octant (1, 1, 1)

   VertexVector m_vertices;   // Vertices [0-8]
                              // Vertex 0 = Corner (0, 0, 0)
                              // Vertex 1 = Corner (1, 0, 0)
                              // Vertex 2 = Corner (0, 0, 1)
                              // Vertex 3 = Corner (1, 0, 1)
                              // Vertex 4 = Corner (0, 1, 0)
                              // Vertex 5 = Corner (1, 1, 0)
                              // Vertex 6 = Corner (0, 1, 1)
                              // Vertex 7 = Corner (1, 1, 1)
                              // Vertex 8 = Center (0.5, 0.5, 0.5)
};

#endif 
