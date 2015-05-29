
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

#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <list>
#include "Common.h"
#include "Util/vectorObj.h"

// Enumerations
enum vertex_color_type
{
   BLACK_VERTEX = 0,
   RED_VERTEX
};

// Vertex class
class Vertex
{
public:
  
   Vertex();
   virtual ~Vertex();

   inline double& x() { return m_coords[0]; }
   inline double& y() { return m_coords[1]; }
   inline double& z() { return m_coords[2]; }

   inline double x() const { return m_coords[0]; }
   inline double y() const { return m_coords[1]; }
   inline double z() const { return m_coords[2]; }
   
   VectorObj&
   getCoords(); 
   void
   setCoords( const VectorObj& coords );

   inline double& grid_x() { return m_grid_coords[0]; }
   inline double& grid_y() { return m_grid_coords[1]; }
   inline double& grid_z() { return m_grid_coords[2]; }

   inline double grid_x() const { return m_grid_coords[0]; }
   inline double grid_y() const { return m_grid_coords[1]; }
   inline double grid_z() const { return m_grid_coords[2]; }

   VectorObj&
   getGridCoords(); 
   void
   setGridCoords( const VectorObj& grid_coords );

   EdgeList&
   getEdgeList();
   void
   setEdgeList( const EdgeList& edge_list );

   TetrahedronVector&
   getTetrahedronVector();
   void
   setTetrahedronVector( const TetrahedronVector& tetrahedron_vector ); 
   
   Edge*
   getEdge( Vertex* end_vertex ) const;

   vertex_color_type
   getColor() const;  
   void
   setColor( vertex_color_type color );

   double
   getF() const;  
   void 
   setF( double f );  
   bool
   getFComputed() const;
   
   bool
   getInSubset() const;
   void
   setInSubset( bool in_subset );

private:

   VectorObj m_coords;                       // 3-D coordinates
   VectorObj m_grid_coords;                  // Grid coordinates: (lowest octant depth - 1) coordinates
   
   EdgeList m_edge_list;                     // List of adjacent edges
   TetrahedronVector m_tetrahedron_vector;   // Adjacent tetrahedra to this vertex
   
   vertex_color_type m_color;                // Vertex color: Black = octant corner, Red = octant center
   double m_f;                               // Value of signed distance function f at this point
   bool m_f_computed;                        // Was the value of f computed?
   bool m_in_subset;                         // Flag to indicate this vertex is in the lattice/octree subset
   
};

#endif 
