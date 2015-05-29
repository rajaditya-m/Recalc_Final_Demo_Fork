
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

#ifndef _OCTREE_H_
#define _OCTREE_H_

#include "Common.h"
#include "Octant.h"
#include "Util/vectorObj.h"

// Octree class
class Octree
{
public:
  
   Octree();
   virtual ~Octree();

   int
   getDepth() const;  
   void 
   setDepth( int depth );
   
   double
   getLeafWidth() const;
   void
   setLeafWidth( double leaf_width );
   
   const VectorObj&
   getOrigin() const;
   void
   setOrigin( const VectorObj& origin );
   
   VertexHashMap&
   getVertices();
   
   EdgeVector&
   getEdges();
      
   Octant*
   getOctant( int depth, int grid_x, int grid_y, int grid_z );
   Octant*
   getOctant( int depth, long hash_key );

   OctantHashMapVector&
   getOctantHashMapVector();
   OctantHashMap&
   getOctantHashMap( int depth );

   Octant*
   createOctant( 
      int depth, int grid_x, int grid_y, int grid_z, 
      std::vector<double>& probe_f, bool probe_f_provided );

   long
   computeOctantHashKey( int depth, int grid_x, int grid_y, int grid_z ) const;
   long
   computeOctantHashKey( int depth, const VectorObj& grid_coords ) const;
   void
   inverseOctantHashKey( int depth, long hash_key, int& grid_x, int& grid_y, int& grid_z ) const;

   long
   computeVertexHashKey( int grid_x, int grid_y, int grid_z ) const;
   long
   computeVertexHashKey( const VectorObj& grid_coords ) const;

   void
   computeVertexGridCoordinates(
      int depth, int octant_grid_x, int octant_grid_y, int octant_grid_z,
      GridCoordVector& vertex_grid_coords );

   bool
   existSmallerAdjacentOctants( int depth, int grid_size, int grid_x, int grid_y, int grid_z );
   bool
   existSmallerAdjacentOctants( int depth, int grid_size, int grid_x, int grid_y, int grid_z, int face );

   Octant*
   getAdjacentOctant( int depth, int grid_size, int grid_x, int grid_y, int grid_z, int face ); 

   Vertex*
   getVertex( int grid_x, int grid_y, int grid_z );
   Vertex*
   getVertex( const VectorObj& grid_coords );
   void
   getEdge( 
      int depth, int grid_size, int grid_x, int grid_y, int grid_z, 
      int edge, Vertex*& p1, Vertex*& p2, Vertex*& center );
   void
   getOctantFace(
      int depth, int grid_size, int grid_x, int grid_y, int grid_z, 
      int face, Vertex*& p1, Vertex*& p2, Vertex*& p3, Vertex*& p4, Vertex*& center );

   void
   getVertexIndicesForFace( int face, VertexIndexVector& vertex_indices );
   void
   getEdgeIndicesForFace( int face, EdgeIndexList& edge_indices );
   
   void
   createEdge( Vertex* v1, Vertex* v2, edge_color_type edge_color ); 
   
   void
   clear();

private:

   int m_depth;                                 // Depth of the octree (0 = root, 1 = 1st level children, etc.)
   double m_leaf_width;                         // Leaf Width - specifies approximate size of tetrahedra. 
   VectorObj m_origin;                          // 3D origin of the octree

   OctantHashMapVector m_hash_map_vector;       // Vector of hash maps for fast lookup of any octant.
                                                // Each hash map represents a depth level of the tree
                                                // (i.e. [0] = {root}, [1] = {child0, child1, ..., child7}, etc.)
                                                // Each hash map is keyed on the grid coordinates of the octant
                                                // within the level.
   
   VertexHashMap m_vertices;                    // Hash map of all vertices in the octree, keyed on the grid
                                                // coordinates.
                                                
   EdgeVector m_edges;                          // Vector of all edges in this octree
         
};

#endif 
