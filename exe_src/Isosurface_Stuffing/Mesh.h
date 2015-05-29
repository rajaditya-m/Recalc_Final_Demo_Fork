
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

#ifndef _MESH_H_
#define _MESH_H_

#include "Common.h"
#include "Octree.h"
#include "Util/Timer.h"
#include <functional>

#define SURFACE_TOLERANCE      1.0E-03        // Tolerance to snap points to surface
#define ALPHA_LONG             0.24999        // Violation parameter for long edges
#define ALPHA_SHORT            0.41189        // Violation parameter for short edges
#define MESH_CUT_Z             15.0           // Cut-away the mesh from the +z or -z direction
#define MESH_CUT               0              // 0 = false, 1 = true
#define MESH_CUT_ZPOSITIVE     1              // 1 = cut in +z direction, 0 = cut in -z direction
#define TETRAHEDRON_TOLERANCE  1.0E-08        // Tolerance for finding point inside tetrahedron
#define MIN_DIHEDRAL_ANGLE     0.162614       // Minimum angle for throwing away bad tetrahedra (in radians)

// Mesh class
class Mesh
{
public:
   
   Mesh();
   virtual ~Mesh();
      
   void 
   generateMesh(std::function<double (const VectorObj &)> phi,
      double octree_depth,
      double octree_leaf_width,
      const VectorObj& octree_origin,
      const VectorObj& starting_octant );
   
   const TetrahedronVector&
   getTetrahedra() const;
   
   void 
   drawMesh();
   void 
   drawWireframeMesh();
   void 
   drawOctreeLeaves();
   void 
   drawOctreeParents();
   void
   drawCutPoints();              
   
   void
   clear();
                                    
public:

   Octree               m_octree;                  // Octree
   TetrahedronVector    m_tetrahedra;              // Tetrahedra
   VertexVector         m_vertices;                // Vertices in the mesh (subset of octree vertices)
   Timer                m_timer;                   // Timer for statistics
   Timer                m_timer2;                  // Timer for statistics

   void 
   drawTetrahedron( const Tetrahedron& tetrahedron );
   
   void 
   computeIntersection(
      std::function<double (const VectorObj& coords)> phi,
      double x1, double y1, double z1,
      double x2, double y2, double z2,
      double& x, double& y, double& z );
   
   void 
   computeCutPoints( 
      std::function<double (const VectorObj& coords)> phi,
      VertexList& vertex_subset );
  
   void 
   warpGrid( VertexList& vertex_subset );
   
   void 
   computeOutputTetrahedra( 
      std::function<double (const VectorObj& coords)> phi,
      const TetrahedronVector& background_grid_tetrahedra );
   
   bool 
   parityRule( Vertex* a, Vertex* b, Vertex* c, Vertex* d );
      
   bool 
   computeSurfaceOctants(std::function<double (const VectorObj& coords)> phi,
      const VectorObj& starting_octant);
   
   void 
   enforceContinuationCondition(std::function<double (const VectorObj& coords)> phi);
   
   void 
   buildOctree();
   
   void 
   enforceWeakBalanceCondition();
   
   void 
   computeBackgroundGrid( 
    std::function<double (const VectorObj& coords)> phi,
      VertexList& vertex_subset,
      TetrahedronVector& background_grid );
   
   void 
   createBackgroundGridTetrahedron(
    std::function<double (const VectorObj& coords)> phi,
      Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4,
      edge_color_type edge_color1, 
      edge_color_type edge_color2, 
      edge_color_type edge_color3, 
      edge_color_type edge_color4, 
      edge_color_type edge_color5, 
      edge_color_type edge_color6,
      tetrahedron_type type,    
      VertexList& vertex_subset, TetrahedronVector& background_grid );

};

#endif 
