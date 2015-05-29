
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

#include <iostream>
#include <math.h>
#include <stdio.h>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <functional>
#include <set>
#include <stack>
#include <cfloat>

#include "Mesh.h"
#include "Util/Util.h"
#include "Vertex.h"
#include "Edge.h"
#include "Tetrahedron.h"


//-----------------------------------------------------
Mesh::Mesh()
{
   // Constructor
   
   m_timer.start();
   m_timer2.start();
}

//-----------------------------------------------------
Mesh::~Mesh()
{
   // Destructor

   clear();
}

//-------------------------------------------------------------------------
void
Mesh::generateMesh( 
    std::function<double (const VectorObj& coords)> phi,
//   double (*phi)(const VectorObj& coords),
   double octree_depth,
   double octree_leaf_width,
   const VectorObj& octree_origin,
   const VectorObj& starting_octant )
{
   // Perform the meshing algorithm

   m_timer.restart();

   // Set octree parameters
   m_octree.setDepth( octree_depth );
   m_octree.setLeafWidth( octree_leaf_width );
   m_octree.setOrigin( octree_origin );   
   
   // Find all leaf octants of the octree that the surface passes through.
   // This step adds leaf octants to the octree.
   std::cout << "Surface ocants calculating.....";
   m_timer2.restart();
   computeSurfaceOctants( phi, starting_octant ); 
   std::cout << m_timer2.getMillisSinceStart() / 1000.0 << "\n";   
         
   // Enforce the continuation condition - add more cubes to ensure only
   // BCC tetrahedra intersect the surface.
   std::cout << "Enforcing continuation condition.....";
   m_timer2.restart();
   enforceContinuationCondition( phi );
   std::cout << m_timer2.getMillisSinceStart() / 1000.0 << "\n";   
   
   // Build the octree (create ancestors of leaf octants)
   std::cout << "Building octree.....";
   m_timer2.restart();
   buildOctree();
   std::cout << m_timer2.getMillisSinceStart() / 1000.0 << "\n";   
   
   // Enforce the weak balance condition - make sure octants only span by at most
   // a factor of 2.
   std::cout << "Enforcing weak balance condition.....";
   m_timer2.restart();
   enforceWeakBalanceCondition(); 
   std::cout << m_timer2.getMillisSinceStart() / 1000.0 << "\n";   
  
   // Convert the balanced octree to a background grid of tetrahedra (slots)
   std::cout << "Computing background grid.....";
   m_timer2.restart();
   VertexList vertex_subset;
   TetrahedronVector background_grid; 
   computeBackgroundGrid( phi, vertex_subset, background_grid );
   std::cout << m_timer2.getMillisSinceStart() / 1000.0 << "\n";   
               
   // For each edge with both end points in the subset, compute
   // an approximate cut point where it intersects the zero-surface.
   std::cout << "Computing cut points.....";
   m_timer2.restart();
   computeCutPoints( phi, vertex_subset );
   std::cout << m_timer2.getMillisSinceStart() / 1000.0 << "\n";   
      
   // Warp the background grid - snap to the surface
   std::cout << "Warping grid.....";
   m_timer2.restart();
   warpGrid( vertex_subset );
   std::cout << m_timer2.getMillisSinceStart() / 1000.0 << "\n";   
      
   // Compute the final output tetrahedra by stuffing the tetrahedra
   // with stencils
   std::cout << "Stuffing tetrahedra.....";
   m_timer2.restart();
   computeOutputTetrahedra( phi, background_grid );
   std::cout << m_timer2.getMillisSinceStart() / 1000.0;
   std::cout << "\nGenerated Mesh in " << m_timer.getMillisSinceStart() / 1000.0 << " seconds.\n";     
   
   // Destroy the background grid
   TetrahedronVectorIterator tet_iter = background_grid.begin();
   for ( ; tet_iter != background_grid.end(); ++tet_iter )
   {
      delete *tet_iter;
   }
   background_grid.clear();  
}

//-------------------------------------------------------------------------
const TetrahedronVector&
Mesh::getTetrahedra() const
{
   return m_tetrahedra;
}

//-------------------------------------------------------------------------
void
Mesh::drawMesh()
{
   // Setup material properties
   GLfloat specular[] = { 0.4, 0.4, 0.4, 1.0 };
   GLfloat shininess[] = { 10.0 };
   glMaterialfv( GL_FRONT, GL_SPECULAR, specular );
   glMaterialfv( GL_FRONT, GL_SHININESS, shininess );
      
   TetrahedronVectorIterator iter = m_tetrahedra.begin();
   for ( ; iter != m_tetrahedra.end(); ++iter )
   {
      Tetrahedron* tetrahedron = *iter;
      
      // If we're cutting the mesh
      if ( MESH_CUT != 0 )
      {
         Vertex* vertex1 = NULL;    
         Vertex* vertex2 = NULL;    
         Vertex* vertex3 = NULL;    
         Vertex* vertex4 = NULL;
         tetrahedron->getVertices( vertex1, vertex2, vertex3, vertex4 );
         
         // Skip this tetrahedron if it has any vertices past the cut-off plane
         if ( MESH_CUT_ZPOSITIVE ) {
            if ( vertex1->z() >= MESH_CUT_Z || vertex2->z() >= MESH_CUT_Z ||
                 vertex3->z() >= MESH_CUT_Z || vertex4->z() >= MESH_CUT_Z )
            {
               continue;
            }
         }
         else {
            if ( vertex1->z() <= MESH_CUT_Z || vertex2->z() <= MESH_CUT_Z ||
                 vertex3->z() <= MESH_CUT_Z || vertex4->z() <= MESH_CUT_Z )
            {
               continue;
            }
         }   
      }

      drawTetrahedron( *tetrahedron );
   }  
}

//-------------------------------------------------------------------------
void
Mesh::drawWireframeMesh()
{
   glColor3d( 0.0, 0.0, 0.0 );
   glLineWidth( 0.5 );

   // Polygon mode to just lines
   glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

   drawMesh();
   
   glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

//-------------------------------------------------------------------------
void
Mesh::drawOctreeLeaves()
{
   // Draw the octree leaf octants

   glColor3d( 0.5, 0.5, 0.5 );
   glLineWidth( 1.0 );
   glBegin( GL_LINES );

   int depth = m_octree.getDepth();
   double grid_size = pow(2.0, depth); 
   
   for ( int x = 0; x < grid_size; ++x )
   {
      for ( int y = 0; y < grid_size; ++y )
      {
         for ( int z = 0; z < grid_size; ++z )
         {
            Octant* octant = m_octree.getOctant( depth, x, y, z );
            if ( octant != NULL )
            {
               VertexVector& vertices = octant->getVertices();
               glVertex3d( vertices[0]->x(), vertices[0]->y(), vertices[0]->z() );
               glVertex3d( vertices[1]->x(), vertices[1]->y(), vertices[1]->z() );

               glVertex3d( vertices[0]->x(), vertices[0]->y(), vertices[0]->z() );
               glVertex3d( vertices[2]->x(), vertices[2]->y(), vertices[2]->z() );

               glVertex3d( vertices[1]->x(), vertices[1]->y(), vertices[1]->z() );
               glVertex3d( vertices[3]->x(), vertices[3]->y(), vertices[3]->z() );

               glVertex3d( vertices[2]->x(), vertices[2]->y(), vertices[2]->z() );
               glVertex3d( vertices[3]->x(), vertices[3]->y(), vertices[3]->z() );

               glVertex3d( vertices[4]->x(), vertices[4]->y(), vertices[4]->z() );
               glVertex3d( vertices[5]->x(), vertices[5]->y(), vertices[5]->z() );

               glVertex3d( vertices[5]->x(), vertices[5]->y(), vertices[5]->z() );
               glVertex3d( vertices[7]->x(), vertices[7]->y(), vertices[7]->z() );

               glVertex3d( vertices[6]->x(), vertices[6]->y(), vertices[6]->z() );
               glVertex3d( vertices[7]->x(), vertices[7]->y(), vertices[7]->z() );

               glVertex3d( vertices[4]->x(), vertices[4]->y(), vertices[4]->z() );
               glVertex3d( vertices[6]->x(), vertices[6]->y(), vertices[6]->z() );

               glVertex3d( vertices[2]->x(), vertices[2]->y(), vertices[2]->z() );
               glVertex3d( vertices[6]->x(), vertices[6]->y(), vertices[6]->z() );

               glVertex3d( vertices[3]->x(), vertices[3]->y(), vertices[3]->z() );
               glVertex3d( vertices[7]->x(), vertices[7]->y(), vertices[7]->z() );

               glVertex3d( vertices[4]->x(), vertices[4]->y(), vertices[4]->z() );
               glVertex3d( vertices[0]->x(), vertices[0]->y(), vertices[0]->z() );

               glVertex3d( vertices[5]->x(), vertices[5]->y(), vertices[5]->z() );
               glVertex3d( vertices[1]->x(), vertices[1]->y(), vertices[1]->z() );
            }
         }
      }
   }
   
   glEnd();
}

//-------------------------------------------------------------------------
void
Mesh::drawOctreeParents()
{
   // Draw the octree parent octants

   int depth = m_octree.getDepth();

   glColor3d( 1.0, 0.0, 1.0 );
   glLineWidth( 1.0 );
   glBegin( GL_LINES );
   
   for ( int current_depth = depth-1; current_depth >= 0; --current_depth )
   {   
      double grid_size = pow( 2.0, current_depth );
   
      for ( int x = 0; x < grid_size; ++x )
      {
         for ( int y = 0; y < grid_size; ++y )
         {
            for ( int z = 0; z < grid_size; ++z )
            {
               Octant* octant = m_octree.getOctant( current_depth, x, y, z );
               if ( octant != NULL )
               {
                  VertexVector& vertices = octant->getVertices();
                  glVertex3d( vertices[0]->x(), vertices[0]->y(), vertices[0]->z() );
                  glVertex3d( vertices[1]->x(), vertices[1]->y(), vertices[1]->z() );

                  glVertex3d( vertices[0]->x(), vertices[0]->y(), vertices[0]->z() );
                  glVertex3d( vertices[2]->x(), vertices[2]->y(), vertices[2]->z() );

                  glVertex3d( vertices[1]->x(), vertices[1]->y(), vertices[1]->z() );
                  glVertex3d( vertices[3]->x(), vertices[3]->y(), vertices[3]->z() );

                  glVertex3d( vertices[2]->x(), vertices[2]->y(), vertices[2]->z() );
                  glVertex3d( vertices[3]->x(), vertices[3]->y(), vertices[3]->z() );

                  glVertex3d( vertices[4]->x(), vertices[4]->y(), vertices[4]->z() );
                  glVertex3d( vertices[5]->x(), vertices[5]->y(), vertices[5]->z() );

                  glVertex3d( vertices[5]->x(), vertices[5]->y(), vertices[5]->z() );
                  glVertex3d( vertices[7]->x(), vertices[7]->y(), vertices[7]->z() );

                  glVertex3d( vertices[6]->x(), vertices[6]->y(), vertices[6]->z() );
                  glVertex3d( vertices[7]->x(), vertices[7]->y(), vertices[7]->z() );

                  glVertex3d( vertices[4]->x(), vertices[4]->y(), vertices[4]->z() );
                  glVertex3d( vertices[6]->x(), vertices[6]->y(), vertices[6]->z() );

                  glVertex3d( vertices[2]->x(), vertices[2]->y(), vertices[2]->z() );
                  glVertex3d( vertices[6]->x(), vertices[6]->y(), vertices[6]->z() );

                  glVertex3d( vertices[3]->x(), vertices[3]->y(), vertices[3]->z() );
                  glVertex3d( vertices[7]->x(), vertices[7]->y(), vertices[7]->z() );

                  glVertex3d( vertices[4]->x(), vertices[4]->y(), vertices[4]->z() );
                  glVertex3d( vertices[0]->x(), vertices[0]->y(), vertices[0]->z() );

                  glVertex3d( vertices[5]->x(), vertices[5]->y(), vertices[5]->z() );
                  glVertex3d( vertices[1]->x(), vertices[1]->y(), vertices[1]->z() );
               }
            }
         }
      }
   }
   
   glEnd();   
}

//-------------------------------------------------------------------------
void
Mesh::drawCutPoints()
{
   EdgeVector& edges = m_octree.getEdges();
   EdgeVectorIterator iter;
   glColor3d(0.0,1.0,0.0);
   glPointSize( 5.0 );
   glBegin(GL_POINTS);
   iter = edges.begin();   
   for ( ; iter != edges.end(); ++iter )
   {
      if ( (*iter)->getCutPoint() != NULL )
         glVertex3d( (*iter)->getCutPoint()->x(), (*iter)->getCutPoint()->y(), (*iter)->getCutPoint()->z() );
   }
   glEnd();
}

//-------------------------------------------------------------------------
void
Mesh::clear()
{
   // Clear the contents of this mesh
   
   m_octree.clear();  
   
   // Destroy the tetrahedra
   TetrahedronVectorIterator tet_iter = m_tetrahedra.begin();
   for ( ; tet_iter != m_tetrahedra.end(); ++tet_iter )
   {
      delete *tet_iter;
   }
   m_tetrahedra.clear();
   
   m_vertices.clear();
}

//-------------------------------------------------------------------------
void 
Mesh::drawTetrahedron( const Tetrahedron& tetrahedron )
{
   GLfloat yellow[] = { 1.0, 1.0, 0.0, 1.0 };
   GLfloat green[] = { 0.0, 1.0, 0.0, 1.0 };

   Vertex* vertex1 = NULL;
   Vertex* vertex2 = NULL;
   Vertex* vertex3 = NULL;
   Vertex* vertex4 = NULL;
   tetrahedron.getVertices( vertex1, vertex2, vertex3, vertex4 );

   glBegin( GL_TRIANGLES );

   // Triangle 1: v1, v2, v3
   if ( vertex1->getF() == 0.0 && vertex2->getF() == 0.0 && vertex3->getF() == 0.0 ) {
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, green );
   }
   else {
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, yellow );
   }

   VectorObj normal = 
      (vertex2->getCoords() - vertex1->getCoords()).cross(
         (vertex3->getCoords() - vertex1->getCoords()) );
   normal.normalize();
   if ( normal.dot( vertex4->getCoords() - vertex1->getCoords() ) > 0.0 ) {
      normal = -normal;
   }
      
   glNormal3d( normal[0], normal[1], normal[2] );
   glVertex3d( vertex1->x(), vertex1->y(), vertex1->z() );
   glNormal3d( normal[0], normal[1], normal[2] );
   glVertex3d( vertex2->x(), vertex2->y(), vertex2->z() );
   glNormal3d( normal[0], normal[1], normal[2] );
   glVertex3d( vertex3->x(), vertex3->y(), vertex3->z() );

   // Triangle 2: v1, v2, v4
   if ( vertex1->getF() == 0.0 && vertex2->getF() == 0.0 && vertex4->getF() == 0.0 ) {
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, green );
   }
   else {
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, yellow );
   }

   normal = 
      (vertex2->getCoords() - vertex1->getCoords()).cross(
         (vertex4->getCoords() - vertex1->getCoords()) );
   normal.normalize();
   if ( normal.dot( vertex3->getCoords() - vertex1->getCoords() ) > 0.0 ) {
      normal = -normal;
   }
   
   glNormal3d( normal[0], normal[1], normal[2] );
   glVertex3d( vertex1->x(), vertex1->y(), vertex1->z() );
   glVertex3d( vertex2->x(), vertex2->y(), vertex2->z() );
   glVertex3d( vertex4->x(), vertex4->y(), vertex4->z() );

   // Triangle 3: v1, v3, v4
   if ( vertex1->getF() == 0.0 && vertex3->getF() == 0.0 && vertex4->getF() == 0.0 ) {
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, green );
   }
   else {
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, yellow );
   }

   normal = 
      (vertex3->getCoords() - vertex1->getCoords()).cross(
         (vertex4->getCoords() - vertex1->getCoords()) );
   normal.normalize();
   if ( normal.dot( vertex2->getCoords() - vertex1->getCoords() ) > 0.0 ) {
      normal = -normal;
   }
      
   glNormal3d( normal[0], normal[1], normal[2] );
   glVertex3d( vertex1->x(), vertex1->y(), vertex1->z() );
   glVertex3d( vertex3->x(), vertex3->y(), vertex3->z() );
   glVertex3d( vertex4->x(), vertex4->y(), vertex4->z() );

   // Triangle 4: v2, v3, v4
   if ( vertex2->getF() == 0.0 && vertex3->getF() == 0.0 && vertex4->getF() == 0.0 ) {
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, green );
   }
   else {
      glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, yellow );
   }

   normal = 
      (vertex3->getCoords() - vertex2->getCoords()).cross(
         (vertex4->getCoords() - vertex2->getCoords()) );
   normal.normalize();
   if ( normal.dot( vertex1->getCoords() - vertex2->getCoords() ) > 0.0 ) {
      normal = -normal;
   }
      
   glNormal3d( normal[0], normal[1], normal[2] );
   glVertex3d( vertex2->x(), vertex2->y(), vertex2->z() );
   glVertex3d( vertex3->x(), vertex3->y(), vertex3->z() );
   glVertex3d( vertex4->x(), vertex4->y(), vertex4->z() );
   
   glEnd();
}

//-------------------------------------------------------------------------
void 
Mesh::computeIntersection(
      std::function<double (const VectorObj& coords)> phi,
   double x1, double y1, double z1,
   double x2, double y2, double z2,
   double& x, double& y, double& z )
{
   // Compute the cut point where the segment between the two given
   // vertices intersects the isosurface.  Use iterative bisection
   // and the cut function phi.  It is assumed the two points are of
   // opposite sign.
         
   double fa = phi( VectorObj(x1, y1, z1) );
   double fp = 0.0;
   
   // Loop enough times
   for ( int i = 0; i < 100; ++i )
   {     
      // Compute the midpoint of the segment
      x = (x1 + x2) * 0.5;
      y = (y1 + y2) * 0.5;
      z = (z1 + z2) * 0.5;

      // If we're within our tolerance, break
      if ( sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) * 0.5 < SURFACE_TOLERANCE ) {
         break;         
      }
            
      // Evaluate function at midpoint
      fp = phi( VectorObj(x, y, z) );
      
      // If we found the intersection, break immediately
      if ( fp == 0.0 ) {
         break;
      }
      
      // Otherwise check if the intersection is in the leftmost interval
      else if ( fa * fp < 0.0 ) 
      {
         x2 = x;
         y2 = y;
         z2 = z;
      }
      
      // Otherwise the intersection is in the rightmost interval
      else 
      {
         x1 = x;
         y1 = y;
         z1 = z;
      
         fa = fp;
      }
   }
}

//-------------------------------------------------------------------------
void 
Mesh::computeCutPoints(
    std::function<double (const VectorObj& coords)> phi,
   VertexList& vertex_subset )
{
   // For each edge that has end points in the subset, compute
   // the cut points that intersect the surface.
   
   // Loop through subset
   Vertex* vertex = NULL;
   Vertex* vertex1 = NULL;
   Vertex* vertex2 = NULL;
   Vertex* cut_point = NULL;
   Edge* edge = NULL;
   double x = 0.0;
   double y = 0.0;
   double z = 0.0;
   VertexListIterator vertex_iter = vertex_subset.begin();
   for ( ; vertex_iter != vertex_subset.end(); ++vertex_iter )
   {
      vertex = *vertex_iter;
      
      // Loop through all adjacent edges
      EdgeList& edge_list = vertex->getEdgeList();
      EdgeListIterator edge_iter = edge_list.begin();
      for ( ; edge_iter != edge_list.end(); ++edge_iter )
      {
         edge = *edge_iter;
      
         // If this edge has already been marked, skip
         if ( edge->getFlag() ) {
            continue;
         }
         
         // Get end points
         vertex1 = edge->getVertex1();
         vertex2 = edge->getVertex2();
         
         // If these are both in the subset and one is inside and one is outside
         // the surface
         if ( vertex1->getInSubset() && vertex2->getInSubset() &&
              vertex1->getF() * vertex2->getF() < 0.0 )
         {
            // Create new cutpoint
            cut_point = new Vertex();
                        
            // Compute intersection
            computeIntersection(
               phi, 
               vertex1->x(), vertex1->y(), vertex1->z(), 
               vertex2->x(), vertex2->y(), vertex2->z(),
               x, y, z );
               
            cut_point->x() = x;
            cut_point->y() = y;
            cut_point->z() = z;
            
            cut_point->setF( 0.0 );
            
            // Add this cut point to the edge (assumes memory management)
            edge->setCutPoint( cut_point );
         }
         
         // Mark the flag to prevent re-processing
         edge->setFlag( true );
      }
   }
}

//-------------------------------------------------------------------------
void 
Mesh::warpGrid( VertexList& vertex_subset )
{
   // Warp the vertices to any cut points that are in violation
   // (i.e. too close).  Essentially, we are snapping the lattice/background grid to 
   // the surface.

   // Loop through subset
   Vertex* vertex = NULL;
   Vertex* vertex1 = NULL;
   Vertex* vertex2 = NULL;
   VectorObj vertex1_coords;
   VectorObj vertex2_coords;
   Vertex* cut_point = NULL;
   Edge* edge = NULL;
   Vertex* candidate_cut_point = NULL;
   double candidate_length = 0.0;
   double edge_length = 0.0;
   double dist_to_cut_point = 0.0;
   VertexListIterator vertex_iter = vertex_subset.begin();
   for ( ; vertex_iter != vertex_subset.end(); ++vertex_iter )
   {
      vertex = *vertex_iter;
      
      // Reset candidate cut point to snap to
      candidate_cut_point = NULL;
      
      // Loop through all adjacent edges
      EdgeList& edge_list = vertex->getEdgeList();
      EdgeListIterator edge_iter = edge_list.begin();
      for ( ; edge_iter != edge_list.end(); ++edge_iter )
      {
         edge = *edge_iter;

         // If there is a cut point on this edge
         if ( edge->getCutPoint() != NULL )
         {             
            // Get the cut point
            cut_point = edge->getCutPoint();

            // Get the end points of the edge and their coordinates
            vertex1 = edge->getVertex1();
            vertex2 = edge->getVertex2();
            vertex1_coords = vertex1->getCoords();
            vertex2_coords = vertex2->getCoords();
            
            // Edge length
            edge_length = (vertex1_coords - vertex2_coords).length();
            
            // Compute distance from vertex to cut point
            dist_to_cut_point = (vertex->getCoords() - cut_point->getCoords()).length();
            
            // If this cut point is in violation of the lattice point
            if ( (edge->getColor() == BLACK_EDGE && 
                  dist_to_cut_point < edge_length * ALPHA_LONG) ||
                 (edge->getColor() == RED_EDGE && 
                  dist_to_cut_point < edge_length * ALPHA_SHORT) )
            {
               // Possibly set new candidate cut point to snap to
               // (i.e. it doesn't exist or this is the closest)
               if ( candidate_cut_point == NULL || dist_to_cut_point < candidate_length )
               {
                  candidate_cut_point = cut_point;
                  candidate_length = dist_to_cut_point;
               }
            }
         }         
      }
      
      // If we have a cut point to snap to
      if ( candidate_cut_point != NULL )
      {
         // Warp the lattice point
         vertex->setCoords( candidate_cut_point->getCoords() );
         
         // Set the cut function value to 0 to indicate it is on the surface
         vertex->setF( 0.0 );
         
         // Discard all cut points on the adjacent edges
         edge_iter = edge_list.begin();
         for ( ; edge_iter != edge_list.end(); ++edge_iter )
         {
            edge = *edge_iter;

            if ( edge->getCutPoint() != NULL )
            {
               delete edge->getCutPoint();
               edge->setCutPoint( NULL );
            }        
         }
      }
   }   
}

//-------------------------------------------------------------------------
void
Mesh::computeOutputTetrahedra( 
    std::function<double (const VectorObj& coords)> phi,
   const TetrahedronVector& background_grid_tetrahedra )
{
   // Compute the final tetrahedra by "stuffing" the background grid 
   // tetrahedra with stencils
   
   // Loop through all background grid tetrahedra
   Vertex* vertex1 = NULL;
   Vertex* vertex2 = NULL;
   Vertex* vertex3 = NULL;
   Vertex* vertex4 = NULL;
   TetrahedronVector quadruple_zero;
   TetrahedronVectorConstIterator tet_iter = background_grid_tetrahedra.begin();
   for ( ; tet_iter != background_grid_tetrahedra.end(); ++tet_iter )
   {
      Tetrahedron* background_grid_tetrahedron = *tet_iter;
      
      if ( background_grid_tetrahedron->getPositiveCount() > 0 ) {
         int x = 0;
      }
      
      // If this tetrahedron has all four vertices on the surface,
      // store it in a list for later processing
      if ( background_grid_tetrahedron->getZeroCount() == 4 )
      {
         quadruple_zero.push_back( background_grid_tetrahedron );
         continue;
      }
            
      // For half, quad, and half pyramids tetrahedra, if it has at least
      // one positive vertex, it becomes an output tetrahedron.
      if ( background_grid_tetrahedron->getType() != BCC_TETRAHEDRON ) 
      {
         if ( background_grid_tetrahedron->getPositiveCount() > 0 ) 
         {
            Tetrahedron* output_tetrahedron1 = new Tetrahedron( *background_grid_tetrahedron );
            output_tetrahedron1->connectVertices( m_vertices );
            m_tetrahedra.push_back( output_tetrahedron1 );    
         }
         continue;
      } 
      
      // Stencil 1: Output full tetrahedron
      if ( background_grid_tetrahedron->getPositiveCount() == 4 )
      {
         Tetrahedron* output_tetrahedron1 = new Tetrahedron( *background_grid_tetrahedron );
         output_tetrahedron1->connectVertices( m_vertices );
         m_tetrahedra.push_back( output_tetrahedron1 );
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 3 && 
                background_grid_tetrahedron->getZeroCount() == 1 )
      {
         // Stencil 2: Output full tetrahedron   
         
         Tetrahedron* output_tetrahedron1 = new Tetrahedron( *background_grid_tetrahedron );
         output_tetrahedron1->connectVertices( m_vertices );
         m_tetrahedra.push_back( output_tetrahedron1 );
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 2 && 
                background_grid_tetrahedron->getZeroCount() == 2 )
      {
         // Stencil 3: Output full tetrahedron   
         
         Tetrahedron* output_tetrahedron1 = new Tetrahedron( *background_grid_tetrahedron );
         output_tetrahedron1->connectVertices( m_vertices );
         m_tetrahedra.push_back( output_tetrahedron1 );
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 1 && 
                background_grid_tetrahedron->getZeroCount() == 3 )
      {
         // Stencil 4: Output full tetrahedron   
         
         Tetrahedron* output_tetrahedron1 = new Tetrahedron( *background_grid_tetrahedron );
         output_tetrahedron1->connectVertices( m_vertices );
         m_tetrahedra.push_back( output_tetrahedron1 );
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 1 && 
                background_grid_tetrahedron->getZeroCount() == 2 &&
                background_grid_tetrahedron->getNegativeCount() == 1 )
      {
         // Stencil 5: Output single truncated tetrahedron   
         
         Vertex* positive_vertex = NULL;
         Vertex* negative_vertex = NULL;
         Vertex* zero_vertex1 = NULL;
         Vertex* zero_vertex2 = NULL;
         
         background_grid_tetrahedron->getVertices( vertex1, vertex2, vertex3, vertex4 );
         
         if ( vertex1->getF() > 0.0 ) positive_vertex = vertex1;
         if ( vertex2->getF() > 0.0 ) positive_vertex = vertex2;
         if ( vertex3->getF() > 0.0 ) positive_vertex = vertex3;
         if ( vertex4->getF() > 0.0 ) positive_vertex = vertex4;

         if ( vertex1->getF() < 0.0 ) negative_vertex = vertex1;
         if ( vertex2->getF() < 0.0 ) negative_vertex = vertex2;
         if ( vertex3->getF() < 0.0 ) negative_vertex = vertex3;
         if ( vertex4->getF() < 0.0 ) negative_vertex = vertex4;
         
         if ( vertex1->getF() == 0.0 ) {
            if ( zero_vertex1 == NULL ) 
               zero_vertex1 = vertex1;
            else
               zero_vertex2 = vertex1;
         }
         if ( vertex2->getF() == 0.0 ) {
            if ( zero_vertex1 == NULL ) 
               zero_vertex1 = vertex2;
            else
               zero_vertex2 = vertex2;
         }
         if ( vertex3->getF() == 0.0 ) {
            if ( zero_vertex1 == NULL ) 
               zero_vertex1 = vertex3;
            else
               zero_vertex2 = vertex3;
         }
         if ( vertex4->getF() == 0.0 ) {
            if ( zero_vertex1 == NULL ) 
               zero_vertex1 = vertex4;
            else
               zero_vertex2 = vertex4;
         }
         
         // Get only cut point
         Vertex* cut_point = NULL;
         Edge* edge = NULL;
         EdgeList& edge_list = positive_vertex->getEdgeList();
         EdgeListIterator edge_iter = edge_list.begin();
         for ( ; edge_iter != edge_list.end(); ++edge_iter )
         {
            edge = *edge_iter;
         
            if ( edge->getVertex1() == positive_vertex && edge->getVertex2() == negative_vertex ||
                 edge->getVertex1() == negative_vertex && edge->getVertex2() == positive_vertex ) 
            {
               cut_point = edge->getCutPoint();
               break;
            }
         }

         Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
         output_tetrahedron1->setVertices( positive_vertex, zero_vertex1, zero_vertex2, cut_point );
         output_tetrahedron1->connectVertices( m_vertices );
         m_tetrahedra.push_back( output_tetrahedron1 );
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 1 && 
                background_grid_tetrahedron->getZeroCount() == 1 &&
                background_grid_tetrahedron->getNegativeCount() == 2 )
      {
         // Stencil 6: Output single truncated tetrahedron   
         
         Vertex* positive_vertex = NULL;
         Vertex* zero_vertex = NULL;
         Vertex* negative_vertex1 = NULL;
         Vertex* negative_vertex2 = NULL;
         
         background_grid_tetrahedron->getVertices( vertex1, vertex2, vertex3, vertex4 );
         
         if ( vertex1->getF() > 0.0 ) positive_vertex = vertex1;
         if ( vertex2->getF() > 0.0 ) positive_vertex = vertex2;
         if ( vertex3->getF() > 0.0 ) positive_vertex = vertex3;
         if ( vertex4->getF() > 0.0 ) positive_vertex = vertex4;

         if ( vertex1->getF() == 0.0 ) zero_vertex = vertex1;
         if ( vertex2->getF() == 0.0 ) zero_vertex = vertex2;
         if ( vertex3->getF() == 0.0 ) zero_vertex = vertex3;
         if ( vertex4->getF() == 0.0 ) zero_vertex = vertex4;
         
         if ( vertex1->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex1;
            else
               negative_vertex2 = vertex1;
         }
         if ( vertex2->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex2;
            else
               negative_vertex2 = vertex2;
         }
         if ( vertex3->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex3;
            else
               negative_vertex2 = vertex3;
         }
         if ( vertex4->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex4;
            else
               negative_vertex2 = vertex4;
         }
                           
         // Get cut point1
         Vertex* cut_point1 = NULL;
         Edge* edge = NULL;
         EdgeList& edge_list1 = negative_vertex1->getEdgeList();
         EdgeListIterator edge_iter = edge_list1.begin();
         for ( ; edge_iter != edge_list1.end(); ++edge_iter )
         {
            edge = *edge_iter;
         
            if ( edge->getVertex1() == positive_vertex && edge->getVertex2() == negative_vertex1 ||
                 edge->getVertex1() == negative_vertex1 && edge->getVertex2() == positive_vertex ) 
            {
               cut_point1 = edge->getCutPoint();
               break;
            }
         }

         // Get cut point2
         Vertex* cut_point2 = NULL;
         EdgeList& edge_list2 = negative_vertex2->getEdgeList();
         edge_iter = edge_list2.begin();
         for ( ; edge_iter != edge_list2.end(); ++edge_iter )
         {
            edge = *edge_iter;
         
            if ( edge->getVertex1() == positive_vertex && edge->getVertex2() == negative_vertex2 ||
                 edge->getVertex1() == negative_vertex2 && edge->getVertex2() == positive_vertex ) 
            {
               cut_point2 = edge->getCutPoint();
               break;
            }
         }
         
         Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
         output_tetrahedron1->setVertices( positive_vertex, cut_point1, cut_point2, zero_vertex );
         output_tetrahedron1->connectVertices( m_vertices );
         m_tetrahedra.push_back( output_tetrahedron1 );
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 1 && 
                background_grid_tetrahedron->getNegativeCount() == 3 )
      {
         // Stencil 7: Output single truncated tetrahedron   
         
         Vertex* positive_vertex = NULL;
         Vertex* negative_vertex1 = NULL;
         Vertex* negative_vertex2 = NULL;
         Vertex* negative_vertex3 = NULL;
         
         background_grid_tetrahedron->getVertices( vertex1, vertex2, vertex3, vertex4 );
         
         if ( vertex1->getF() > 0.0 ) positive_vertex = vertex1;
         if ( vertex2->getF() > 0.0 ) positive_vertex = vertex2;
         if ( vertex3->getF() > 0.0 ) positive_vertex = vertex3;
         if ( vertex4->getF() > 0.0 ) positive_vertex = vertex4;

         if ( vertex1->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex1;
            else if ( negative_vertex2 == NULL )
               negative_vertex2 = vertex1;
            else
               negative_vertex3 = vertex1;
         }
         if ( vertex2->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex2;
            else if ( negative_vertex2 == NULL )
               negative_vertex2 = vertex2;
            else
               negative_vertex3 = vertex2;
         }
         if ( vertex3->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex3;
            else if ( negative_vertex2 == NULL )
               negative_vertex2 = vertex3;
            else
               negative_vertex3 = vertex3;
         }
         if ( vertex4->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex4;
            else if ( negative_vertex2 == NULL )
               negative_vertex2 = vertex4;
            else
               negative_vertex3 = vertex4;
         }
                                    
         // Get cut points
         Vertex* cut_point1 = NULL;
         Vertex* cut_point2 = NULL;
         Vertex* cut_point3 = NULL;
         Edge* edge = NULL;
         EdgeList& edge_list = positive_vertex->getEdgeList();
         EdgeListIterator edge_iter = edge_list.begin();
         for ( ; edge_iter != edge_list.end(); ++edge_iter )
         {
            edge = *edge_iter;
         
            if ( edge->getVertex1() == positive_vertex && edge->getVertex2() == negative_vertex1 ||
                 edge->getVertex1() == negative_vertex1 && edge->getVertex2() == positive_vertex ) 
            {
               cut_point1 = edge->getCutPoint();
            }
            else if ( edge->getVertex1() == positive_vertex && edge->getVertex2() == negative_vertex2 ||
                 edge->getVertex1() == negative_vertex2 && edge->getVertex2() == positive_vertex ) 
            {
               cut_point2 = edge->getCutPoint();
            }
            else if ( edge->getVertex1() == positive_vertex && edge->getVertex2() == negative_vertex3 ||
                 edge->getVertex1() == negative_vertex3 && edge->getVertex2() == positive_vertex ) 
            {
               cut_point3 = edge->getCutPoint();
            }
         }
                  
         Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
         output_tetrahedron1->setVertices( positive_vertex, cut_point1, cut_point2, cut_point3 );
         output_tetrahedron1->connectVertices( m_vertices );
         m_tetrahedra.push_back( output_tetrahedron1 );
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 2 && 
                background_grid_tetrahedron->getNegativeCount() == 1 &&
                background_grid_tetrahedron->getZeroCount() == 1 )
      {
         // Stencil 8 and 10: Output two tetrahedron   
         
         Vertex* zero_vertex = NULL;
         Vertex* negative_vertex = NULL;
         Vertex* positive_vertex1 = NULL;
         Vertex* positive_vertex2 = NULL;
         
         background_grid_tetrahedron->getVertices( vertex1, vertex2, vertex3, vertex4 );
         
         if ( vertex1->getF() == 0.0 ) zero_vertex = vertex1;
         if ( vertex2->getF() == 0.0 ) zero_vertex = vertex2;
         if ( vertex3->getF() == 0.0 ) zero_vertex = vertex3;
         if ( vertex4->getF() == 0.0 ) zero_vertex = vertex4;

         if ( vertex1->getF() < 0.0 ) negative_vertex = vertex1;
         if ( vertex2->getF() < 0.0 ) negative_vertex = vertex2;
         if ( vertex3->getF() < 0.0 ) negative_vertex = vertex3;
         if ( vertex4->getF() < 0.0 ) negative_vertex = vertex4;

         if ( vertex1->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex1;
            else
               positive_vertex2 = vertex1;
         }
         if ( vertex2->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex2;
            else
               positive_vertex2 = vertex2;
         }
         if ( vertex3->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex3;
            else
               positive_vertex2 = vertex3;
         }
         if ( vertex4->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex4;
            else
               positive_vertex2 = vertex4;
         }
         
         // Get + to + edge
         Edge* plus_to_plus_edge = NULL;
         EdgeList& edge_list = positive_vertex1->getEdgeList();
         EdgeListIterator edge_iter = edge_list.begin();
         for ( ; edge_iter != edge_list.end(); ++edge_iter )
         {
            plus_to_plus_edge = *edge_iter;
            
            if ( plus_to_plus_edge->getVertex1() == positive_vertex1 && 
                 plus_to_plus_edge->getVertex2() == positive_vertex2 ||
                 plus_to_plus_edge->getVertex1() == positive_vertex2 && 
                 plus_to_plus_edge->getVertex2() == positive_vertex1 ) 
            {
               break;
            }
         }         
                                    
         // Determine if + to + edge is red (i.e. stencil 8)
         if ( plus_to_plus_edge->getColor() == RED_EDGE )
         {
            // Get black edge with cut point
            Vertex* end_vertex = NULL;
            Edge* black_edge = background_grid_tetrahedron->getEdge( negative_vertex, BLACK_EDGE, 1.0, end_vertex );         
         
            // Mark positive vertex 1 and store cut point 1
            positive_vertex1 = end_vertex;
            Vertex* cut_point1 = black_edge->getCutPoint();
            
            // Get red edge with cut point
            Edge* red_edge = background_grid_tetrahedron->getEdge( negative_vertex, RED_EDGE, 1.0, end_vertex );                     

            // Mark positive vertex 2 and store cut point 2
            positive_vertex2 = end_vertex;
            Vertex* cut_point2 = red_edge->getCutPoint();
            
            // Output both tetrahedra
            Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
            Tetrahedron* output_tetrahedron2 = new Tetrahedron(); 
            output_tetrahedron1->setVertices( positive_vertex1, positive_vertex2, cut_point1, zero_vertex );
            output_tetrahedron2->setVertices( cut_point1, cut_point2, positive_vertex2, zero_vertex );
            output_tetrahedron1->connectVertices( m_vertices );
            output_tetrahedron2->connectVertices( m_vertices );
            m_tetrahedra.push_back( output_tetrahedron1 );            
            m_tetrahedra.push_back( output_tetrahedron2 );            
         }
         else
         {
            // Stencil 10: + to + edge is black
            
            // Arbitrarily pick a, b, c, and d vertices for parity rule
            Vertex* end_vertex = NULL;    
            Edge* edge1 = background_grid_tetrahedron->getEdge( positive_vertex1, RED_EDGE, -1.0, end_vertex );
            Vertex* vertex_a = positive_vertex1;
            Vertex* vertex_d = edge1->getCutPoint();
            Edge* edge2 = background_grid_tetrahedron->getEdge( positive_vertex2, RED_EDGE, -1.0, end_vertex );
            Vertex* vertex_b = positive_vertex2;
            Vertex* vertex_c = edge2->getCutPoint();
            
            // Apply parity rule - returns whether ac is our diagonal
            if ( parityRule( vertex_a, vertex_b, vertex_c, vertex_d ) )
            {
               // Output both tetrahedra
               Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
               Tetrahedron* output_tetrahedron2 = new Tetrahedron(); 
               output_tetrahedron1->setVertices( zero_vertex, vertex_d, vertex_a, vertex_c );
               output_tetrahedron2->setVertices( zero_vertex, vertex_a, vertex_c, vertex_b );
               output_tetrahedron1->connectVertices( m_vertices );
               output_tetrahedron2->connectVertices( m_vertices );
               m_tetrahedra.push_back( output_tetrahedron1 );            
               m_tetrahedra.push_back( output_tetrahedron2 );            
            }
            else
            {
               // Output both tetrahedra
               Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
               Tetrahedron* output_tetrahedron2 = new Tetrahedron(); 
               output_tetrahedron1->setVertices( zero_vertex, vertex_d, vertex_c, vertex_b );
               output_tetrahedron2->setVertices( zero_vertex, vertex_d, vertex_a, vertex_b );
               output_tetrahedron1->connectVertices( m_vertices );
               output_tetrahedron2->connectVertices( m_vertices );
               m_tetrahedra.push_back( output_tetrahedron1 );            
               m_tetrahedra.push_back( output_tetrahedron2 );
            } 
         }
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 2 && 
                background_grid_tetrahedron->getNegativeCount() == 2 )
      {
         // Stencil 9 and 12: Output three tetrahedron
         
         Vertex* negative_vertex1 = NULL;
         Vertex* negative_vertex2 = NULL;
         Vertex* positive_vertex1 = NULL;
         Vertex* positive_vertex2 = NULL;
         
         background_grid_tetrahedron->getVertices( vertex1, vertex2, vertex3, vertex4 );

         if ( vertex1->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex1;
            else
               negative_vertex2 = vertex1;
         }
         if ( vertex2->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex2;
            else
               negative_vertex2 = vertex2;
         }
         if ( vertex3->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex3;
            else
               negative_vertex2 = vertex3;
         }
         if ( vertex4->getF() < 0.0 ) {
            if ( negative_vertex1 == NULL ) 
               negative_vertex1 = vertex4;
            else
               negative_vertex2 = vertex4;
         }
                           
         if ( vertex1->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex1;
            else
               positive_vertex2 = vertex1;
         }
         if ( vertex2->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex2;
            else
               positive_vertex2 = vertex2;
         }
         if ( vertex3->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex3;
            else
               positive_vertex2 = vertex3;
         }
         if ( vertex4->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex4;
            else
               positive_vertex2 = vertex4;
         }
         
         // If a red edge exists between the positive vertices (i.e Stencil 9)
         Vertex* end_vertex = NULL;
         Edge* positive_edge = background_grid_tetrahedron->getEdge( positive_vertex1, RED_EDGE, 1.0, end_vertex );
         if ( positive_edge != NULL )
         {
            // This case is symmetrical, but does not require the parity 
            // rule for disambiguation.  Therefore, we arbitrarily choose to assign
            // vertices to letter variables, just as long as the symmetry is preserved.
            // 
            // Assign a and d to the cut points on the + to - red edges
            // Assign b and c to the cut points on the + to - black edges
            // Assign e to the positive vertex adjacent to a and b
            // Assign f to the positive vertex adjacent to c and d
            
            Vertex* vertex_e = positive_vertex1;
            Vertex* vertex_f = positive_vertex2;
            
            Edge* plus_to_minus_edge = 
               background_grid_tetrahedron->getEdge( vertex_e, RED_EDGE, -1.0, end_vertex );
            Vertex* vertex_a = plus_to_minus_edge->getCutPoint();
            
            plus_to_minus_edge = 
               background_grid_tetrahedron->getEdge( vertex_e, BLACK_EDGE, -1.0, end_vertex );
            Vertex* vertex_b = plus_to_minus_edge->getCutPoint();

            plus_to_minus_edge = 
               background_grid_tetrahedron->getEdge( vertex_f, RED_EDGE, -1.0, end_vertex );
            Vertex* vertex_d = plus_to_minus_edge->getCutPoint();

            plus_to_minus_edge = 
               background_grid_tetrahedron->getEdge( vertex_f, BLACK_EDGE, -1.0, end_vertex );
            Vertex* vertex_c = plus_to_minus_edge->getCutPoint();
            
            // Output three tetrahedra
            Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
            Tetrahedron* output_tetrahedron2 = new Tetrahedron(); 
            Tetrahedron* output_tetrahedron3 = new Tetrahedron(); 
            output_tetrahedron1->setVertices( vertex_b, vertex_d, vertex_c, vertex_f );
            output_tetrahedron2->setVertices( vertex_e, vertex_b, vertex_a, vertex_c );
            output_tetrahedron3->setVertices( vertex_e, vertex_b, vertex_c, vertex_f );
            output_tetrahedron1->connectVertices( m_vertices );
            output_tetrahedron2->connectVertices( m_vertices );
            output_tetrahedron3->connectVertices( m_vertices );
            m_tetrahedra.push_back( output_tetrahedron1 );            
            m_tetrahedra.push_back( output_tetrahedron2 );            
            m_tetrahedra.push_back( output_tetrahedron3 );              
         }
         else
         {
            // Stencil 12
            
            // Arbitrarily pick a, b, c, and d vertices for parity rule
            Vertex* vertex_a = positive_vertex1;
            Vertex* vertex_b = positive_vertex2;
            
            // We are choosing an arbitrary edge, as long as we stay symmetrical.
            Edge* edge1 = vertex_a->getEdge( negative_vertex1 );
            Vertex* vertex_d = edge1->getCutPoint();
            Edge* edge2 = vertex_b->getEdge( negative_vertex1 );
            Vertex* vertex_c = edge2->getCutPoint();
            
            // Get the two cut points (which are on the other side opposite of d and c)
            edge1 = vertex_a->getEdge( negative_vertex2 );
            Vertex* cut_point1 = edge1->getCutPoint();
            edge2 = vertex_b->getEdge( negative_vertex2 );
            Vertex* cut_point2 = edge2->getCutPoint();
            
            // Apply parity rule - returns whether ac is our diagonal
            if ( parityRule( vertex_a, vertex_b, vertex_c, vertex_d ) )
            {
               // Output three tetrahedra
               Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
               Tetrahedron* output_tetrahedron2 = new Tetrahedron(); 
               Tetrahedron* output_tetrahedron3 = new Tetrahedron(); 
               output_tetrahedron1->setVertices( cut_point1, vertex_d, vertex_c, vertex_a );
               output_tetrahedron2->setVertices( cut_point1, cut_point2, vertex_c, vertex_b );
               output_tetrahedron3->setVertices( cut_point1, vertex_c, vertex_a, vertex_b );
               output_tetrahedron1->connectVertices( m_vertices );
               output_tetrahedron2->connectVertices( m_vertices );
               output_tetrahedron3->connectVertices( m_vertices );
               m_tetrahedra.push_back( output_tetrahedron1 );            
               m_tetrahedra.push_back( output_tetrahedron2 );            
               m_tetrahedra.push_back( output_tetrahedron3 );            
            }
            else
            {
               // Output three tetrahedra
               Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
               Tetrahedron* output_tetrahedron2 = new Tetrahedron(); 
               Tetrahedron* output_tetrahedron3 = new Tetrahedron(); 
               output_tetrahedron1->setVertices( vertex_d, vertex_c, vertex_b, cut_point2 );
               output_tetrahedron2->setVertices( vertex_d, vertex_a, vertex_b, cut_point2 );
               output_tetrahedron3->setVertices( vertex_d, vertex_a, cut_point1, cut_point2 );
               output_tetrahedron1->connectVertices( m_vertices );
               output_tetrahedron2->connectVertices( m_vertices );
               output_tetrahedron3->connectVertices( m_vertices );
               m_tetrahedra.push_back( output_tetrahedron1 );            
               m_tetrahedra.push_back( output_tetrahedron2 );            
               m_tetrahedra.push_back( output_tetrahedron3 );            
            }            
         }
      }
      else if ( background_grid_tetrahedron->getPositiveCount() == 3 && 
                background_grid_tetrahedron->getNegativeCount() == 1 )
      {
         // Stencil 11: Output three tetrahedron
         
         Vertex* negative_vertex = NULL;
         Vertex* positive_vertex1 = NULL;
         Vertex* positive_vertex2 = NULL;
         Vertex* positive_vertex3 = NULL;
         
         background_grid_tetrahedron->getVertices( vertex1, vertex2, vertex3, vertex4 );
         
         if ( vertex1->getF() < 0.0 ) negative_vertex = vertex1;
         if ( vertex2->getF() < 0.0 ) negative_vertex = vertex2;
         if ( vertex3->getF() < 0.0 ) negative_vertex = vertex3;
         if ( vertex4->getF() < 0.0 ) negative_vertex = vertex4;

         if ( vertex1->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex1;
            else if ( positive_vertex2 == NULL )
               positive_vertex2 = vertex1;
            else
               positive_vertex3 = vertex1;
         }
         if ( vertex2->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex2;
            else if ( positive_vertex2 == NULL )
               positive_vertex2 = vertex2;
            else
               positive_vertex3 = vertex2;
         }
         if ( vertex3->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex3;
            else if ( positive_vertex2 == NULL )
               positive_vertex2 = vertex3;
            else
               positive_vertex3 = vertex3;
         }
         if ( vertex4->getF() > 0.0 ) {
            if ( positive_vertex1 == NULL ) 
               positive_vertex1 = vertex4;
            else if ( positive_vertex2 == NULL )
               positive_vertex2 = vertex4;
            else
               positive_vertex3 = vertex4;
         }
                                             
         // Get the cut point that is not part of the quadrilateral (i.e. on the black edge)
         // Also get the lone positive vertex that is part of this edge
         Vertex* end_vertex = NULL;    
         Edge* neg_black_edge = background_grid_tetrahedron->getEdge( negative_vertex, BLACK_EDGE, 1.0, end_vertex ); 
         Vertex* cut_point = neg_black_edge->getCutPoint();
         Vertex* excluded_positive = end_vertex;
                  
         // Determine which two positive vertices are part of the quadrilateral
         Vertex* quad_positive1 = NULL;
         Vertex* quad_positive2 = NULL;
         if ( end_vertex == positive_vertex1 ) {
            quad_positive1 = positive_vertex2;
            quad_positive2 = positive_vertex3;
         }
         else if ( end_vertex == positive_vertex2 ) {
            quad_positive1 = positive_vertex1;
            quad_positive2 = positive_vertex3;
         }
         else {
            quad_positive1 = positive_vertex1;
            quad_positive2 = positive_vertex2;
         }
                  
         // Arbitrarily pick a, b, c, and d vertices for parity rule
         Edge* edge1 = negative_vertex->getEdge( quad_positive1 );
         Vertex* vertex_a = quad_positive1;
         Vertex* vertex_d = edge1->getCutPoint();
         Edge* edge2 = negative_vertex->getEdge( quad_positive2 );
         Vertex* vertex_b = quad_positive2;
         Vertex* vertex_c = edge2->getCutPoint();
         
         // Apply parity rule - returns whether ac is our diagonal
         if ( parityRule( vertex_a, vertex_b, vertex_c, vertex_d ) )
         {
            // Output three tetrahedra
            Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
            Tetrahedron* output_tetrahedron2 = new Tetrahedron(); 
            Tetrahedron* output_tetrahedron3 = new Tetrahedron(); 
            output_tetrahedron1->setVertices( cut_point, vertex_d, vertex_a, vertex_c );
            output_tetrahedron2->setVertices( cut_point, vertex_c, vertex_a, vertex_b );
            output_tetrahedron3->setVertices( excluded_positive, cut_point, vertex_a, vertex_b );
            output_tetrahedron1->connectVertices( m_vertices );
            output_tetrahedron2->connectVertices( m_vertices );
            output_tetrahedron3->connectVertices( m_vertices );
            m_tetrahedra.push_back( output_tetrahedron1 );            
            m_tetrahedra.push_back( output_tetrahedron2 );            
            m_tetrahedra.push_back( output_tetrahedron3 );            
         }
         else
         {
            // Output three tetrahedra
            Tetrahedron* output_tetrahedron1 = new Tetrahedron(); 
            Tetrahedron* output_tetrahedron2 = new Tetrahedron(); 
            Tetrahedron* output_tetrahedron3 = new Tetrahedron(); 
            output_tetrahedron1->setVertices( cut_point, vertex_d, vertex_c, vertex_b );
            output_tetrahedron2->setVertices( cut_point, vertex_a, vertex_b, vertex_d );
            output_tetrahedron3->setVertices( excluded_positive, cut_point, vertex_a, vertex_b );
            output_tetrahedron1->connectVertices( m_vertices );
            output_tetrahedron2->connectVertices( m_vertices );
            output_tetrahedron3->connectVertices( m_vertices );
            m_tetrahedra.push_back( output_tetrahedron1 );            
            m_tetrahedra.push_back( output_tetrahedron2 );            
            m_tetrahedra.push_back( output_tetrahedron3 );            
         }          
      }      
   }
   
   // Special processing for all quadruple zero tetrahedra
   tet_iter = quadruple_zero.begin();
   for ( ; tet_iter != quadruple_zero.end(); ++tet_iter )
   {
      Tetrahedron* tetrahedron = *tet_iter;
      
      // Get vertices
      Vertex* vertex1 = NULL;
      Vertex* vertex2 = NULL;
      Vertex* vertex3 = NULL;
      Vertex* vertex4 = NULL;
      tetrahedron->getVertices( vertex1, vertex2, vertex3, vertex4 );
      
      // Skip if this tetrahedron has poor dihedral angles
      if ( tetrahedron->isWellConditioned( MIN_DIHEDRAL_ANGLE, TETRAHEDRON_TOLERANCE ) == false ) {
         continue;
      }
                  
      // Skip if all four faces do not adjoin output tetrahedra
      // (Not including the quadruple zero kind)
      
      // Face: V1, V2, V3
      bool adjacent_tet_found1 = false;
      bool adjacent_tet_found2 = false;
      bool adjacent_tet_found3 = false;
      bool adjacent_tet_found4 = false;
      TetrahedronVector& tets1 = vertex1->getTetrahedronVector();
      TetrahedronVectorIterator vt_iter = tets1.begin();
      for ( ; vt_iter != tets1.end(); ++vt_iter )
      {
         Tetrahedron* other_tet = *vt_iter;
         
         // If we found the adjacent tet to this face
         if ( other_tet->containsVertex( vertex1 ) &&
              other_tet->containsVertex( vertex2 ) &&
              other_tet->containsVertex( vertex3 ) ) 
         {
            // Make sure this is not the quadruple kind
            if ( other_tet->getZeroCount() != 4 ) {
               adjacent_tet_found1 = true;
            }        
            break;      
         }
      }

      if ( adjacent_tet_found1 == false )
      {
         // Face: V1, V2, V4
         vt_iter = tets1.begin();
         for ( ; vt_iter != tets1.end(); ++vt_iter )
         {
            Tetrahedron* other_tet = *vt_iter;
            
            // If we found the adjacent tet to this face
            if ( other_tet->containsVertex( vertex1 ) &&
                 other_tet->containsVertex( vertex2 ) &&
                 other_tet->containsVertex( vertex4 ) ) 
            {
               // Make sure this is not the quadruple kind
               if ( other_tet->getZeroCount() != 4 ) {
                  adjacent_tet_found2 = true;
               }        
               break;      
            }
         }
         
         if ( adjacent_tet_found2 == false )
         {
            // Face: V1, V3, V4
            vt_iter = tets1.begin();
            for ( ; vt_iter != tets1.end(); ++vt_iter )
            {
               Tetrahedron* other_tet = *vt_iter;
               
               // If we found the adjacent tet to this face
               if ( other_tet->containsVertex( vertex1 ) &&
                    other_tet->containsVertex( vertex3 ) &&
                    other_tet->containsVertex( vertex4 ) ) 
               {
                  // Make sure this is not the quadruple kind
                  if ( other_tet->getZeroCount() != 4 ) {
                     adjacent_tet_found3 = true;
                  }        
                  break;      
               }
            }

            if ( adjacent_tet_found3 == false )
            {
               // Face: V2, V3, V3
               TetrahedronVector& tets2 = vertex2->getTetrahedronVector();
               vt_iter = tets2.begin();
               for ( ; vt_iter != tets2.end(); ++vt_iter )
               {
                  Tetrahedron* other_tet = *vt_iter;
                  
                  // If we found the adjacent tet to this face
                  if ( other_tet->containsVertex( vertex2 ) &&
                       other_tet->containsVertex( vertex3 ) &&
                       other_tet->containsVertex( vertex4 ) ) 
                  {
                     // Make sure this is not the quadruple kind
                     if ( other_tet->getZeroCount() != 4 ) {
                        adjacent_tet_found4 = true;
                     }        
                     break;      
                  }
               }
            }
         }
      }
      
      // If we did not find any adjacent tets, skip
      if ( adjacent_tet_found1 == false && 
           adjacent_tet_found2 == false && 
           adjacent_tet_found3 == false && 
           adjacent_tet_found4 == false )
      {
         continue;
      }
      
      // Compute tetrahedron center point
      VectorObj center = 
         (vertex1->getCoords() + 
          vertex2->getCoords() + 
          vertex3->getCoords() + 
          vertex4->getCoords()) / 4.0;
          
      // Skip if phi value is outside of mesh
      double f_value = phi( center );
      if ( f_value < 0.0 ) {
         continue;
      }
      
      // Add this quadruple zero tet to the mesh
      Tetrahedron* output_tetrahedron = new Tetrahedron( *tetrahedron );
      output_tetrahedron->connectVertices( m_vertices );
      m_tetrahedra.push_back( output_tetrahedron );          
   } 
}

//-------------------------------------------------------------------------
bool 
Mesh::parityRule( Vertex* a, Vertex* b, Vertex* c, Vertex* d )
{
   // The parity rule to break symmetry when attempting to split
   // a truncated triangle that is a quadrilateral.
   // This rule only applies when the black edge of the triangle
   // has not been split.
   // A and B are the end points of the black edge, and C and D
   // are the cut points of the red edges such that the diagonals
   // are ac and bd.  
   // Returns true if ac is the diagonal, false for bd.
   
   // Determine if a has an even or odd number of coordinates greater
   // than c.
   bool even = true;
   if ( a->x() > c->x() ) {
      even = ! even;
   }
   if ( a->y() > c->y() ) {
      even = ! even;
   }
   if ( a->z() > c->z() ) {
      even = ! even;
   }
      
   // Opposite rules for black vs. red vertices
   if ( a->getColor() == BLACK_VERTEX )
   {
      // Odd = ac diagonal, Even = bd diagonal
      return ! even;
   }
   else
   {
      // Even = ac diagonal, Odd = bd diagonal
      return even;
   }
}

//-------------------------------------------------------------------------
bool 
Mesh::computeSurfaceOctants(
    std::function<double (const VectorObj& coord)> phi,
   const VectorObj& starting_octant)
{
   // Compute the leaf octants of the octree that the surface passes through.
   // Returns true if successful, false if an error occurred 
   // (i.e. initial surface octant not correct)
   
   // Get the octree parameters
   int depth = m_octree.getDepth();
   double leaf_width = m_octree.getLeafWidth();
   VectorObj origin = m_octree.getOrigin();
         
   // Determine dimension of the leaf grid, and its spatial size
   double grid_size = pow( 2.0, depth );
   double size = leaf_width * grid_size;

   // Find the grid coordinates of the initial octant
   int grid_x = (int) (((starting_octant[0] - origin[0]) + (size / 2.0)) / leaf_width); 
   int grid_y = (int) (((starting_octant[1] - origin[1]) + (size / 2.0)) / leaf_width); 
   int grid_z = (int) (((starting_octant[2] - origin[2]) + (size / 2.0)) / leaf_width); 
   
   // Make sure this is in the extent of the octree
   bool not_in_grid = false;
   if ( grid_x < 0 || grid_x >= grid_size || 
        grid_y < 0 || grid_y >= grid_size ||
        grid_z < 0 || grid_z >= grid_size )
   {
      not_in_grid = true;
   }

   // Determine 3D coordinate range of probe points
   double x1 = origin[0] - (size / 2.0) + (grid_x * leaf_width);
   double x2 = origin[0] - (size / 2.0) + ((grid_x+1) * leaf_width);
   double y1 = origin[1] - (size / 2.0) + (grid_y * leaf_width);
   double y2 = origin[1] - (size / 2.0) + ((grid_y+1) * leaf_width);
   double z1 = origin[2] - (size / 2.0) + (grid_z * leaf_width);
   double z2 = origin[2] - (size / 2.0) + ((grid_z+1) * leaf_width);
   
   // Determine cut function values at probe points
   std::vector<double> probe_f;
   probe_f.resize( 9 );
   probe_f[0] = phi( VectorObj(x1, y1, z1) ); 
   probe_f[1] = phi( VectorObj(x2, y1, z1) ); 
   probe_f[2] = phi( VectorObj(x1, y1, z2) ); 
   probe_f[3] = phi( VectorObj(x2, y1, z2) ); 
   probe_f[4] = phi( VectorObj(x1, y2, z1) ); 
   probe_f[5] = phi( VectorObj(x2, y2, z1) ); 
   probe_f[6] = phi( VectorObj(x1, y2, z2) ); 
   probe_f[7] = phi( VectorObj(x2, y2, z2) ); 
   probe_f[8] = phi( VectorObj((x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0) ); 
   
   // Determine if it is on the surface:
   // Does it contain a non-negative and non-positive vertex?  
   // Zero counts as both.
   bool non_negative_found = false;
   bool non_positive_found = false;
   std::vector<double>::iterator probe_iter = probe_f.begin();
   for ( ; probe_iter != probe_f.end(); ++probe_iter )
   {
      if ( *probe_iter >= 0.0 ) {
         non_negative_found = true;
      }
      
      if ( *probe_iter <= 0.0 ) {
         non_positive_found = true;
      }
      
      if ( non_negative_found && non_positive_found ) {
         break;
      }
   }
   
   // If it is not on the surface, search for an initial surface octant
   if ( not_in_grid || non_negative_found == false || non_positive_found == false )
   {
   	std::cout << "\nWarning: Initial surface octant is not on the surface!\n";
   	std::cout << "Warning: Searching manually for surface octant.....";
    
      // Find an initial surface octant
      bool surface_octant_found = false;
      for ( grid_x = 0; grid_x < grid_size; ++grid_x ) 
      {
         for ( grid_y = 0; grid_y < grid_size; ++grid_y ) 
         {
            for ( grid_z = 0; grid_z < grid_size; ++grid_z ) 
            {
               // Determine 3D coordinate range of probe points
               x1 = origin[0] - (size / 2.0) + (grid_x * leaf_width);
               x2 = origin[0] - (size / 2.0) + ((grid_x+1) * leaf_width);
               y1 = origin[1] - (size / 2.0) + (grid_y * leaf_width);
               y2 = origin[1] - (size / 2.0) + ((grid_y+1) * leaf_width);
               z1 = origin[2] - (size / 2.0) + (grid_z * leaf_width);
               z2 = origin[2] - (size / 2.0) + ((grid_z+1) * leaf_width);
               
               // Determine cut function values at probe points
               probe_f[0] = phi( VectorObj(x1, y1, z1) ); 
               probe_f[1] = phi( VectorObj(x2, y1, z1) ); 
               probe_f[2] = phi( VectorObj(x1, y1, z2) ); 
               probe_f[3] = phi( VectorObj(x2, y1, z2) ); 
               probe_f[4] = phi( VectorObj(x1, y2, z1) ); 
               probe_f[5] = phi( VectorObj(x2, y2, z1) ); 
               probe_f[6] = phi( VectorObj(x1, y2, z2) ); 
               probe_f[7] = phi( VectorObj(x2, y2, z2) ); 
               probe_f[8] = phi( VectorObj((x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0) ); 
               
               // Determine if it is on the surface:
               // Does it contain a non-negative and non-positive vertex?  
               // Zero counts as both.
               non_negative_found = false;
               non_positive_found = false;
               probe_iter = probe_f.begin();
               for ( ; probe_iter != probe_f.end(); ++probe_iter )
               {
                  if ( *probe_iter >= 0.0 ) {
                     non_negative_found = true;
                  }
                  
                  if ( *probe_iter <= 0.0 ) {
                     non_positive_found = true;
                  }
                  
                  if ( non_negative_found && non_positive_found ) {
                     break;
                  }
               }
               
               // If it is on the surface, break
               if ( non_negative_found && non_positive_found ) 
               {
                  surface_octant_found = true;
                  break;
               }
            }
            if ( surface_octant_found ) break;
         }
         if ( surface_octant_found ) break;
      }
      
      // Return error if we didn't find the surface octant
      if ( surface_octant_found == false ) {
         return false;
      }
      
      std::cout << "Surface Octant Found!\n";      
   }
   
   // Create a set of hash keys (coordinates) of visited octants for depth first search.
   std::set<long> visited;
   
   // Create a stack of hash keys (coordinates) of pending cubes for depth first search.
   // Initialize with the starting octant.
   std::stack<long> pending;
   pending.push( m_octree.computeOctantHashKey( depth, grid_x, grid_y, grid_z ) );

   // Perform depth first search over all surface cubes
   while ( pending.empty() == false )
   {
      // Pop the current octant
      long hash_key = pending.top();
      pending.pop();
   
      // Mark that this octant has been visited
      visited.insert( hash_key );
      
      // Get grid coordinates
      m_octree.inverseOctantHashKey( depth, hash_key, grid_x, grid_y, grid_z );   

      // Determine 3D coordinate range of probe points for current octant
      x1 = origin[0] - (size / 2.0) + (grid_x * leaf_width);
      x2 = origin[0] - (size / 2.0) + ((grid_x+1) * leaf_width);
      y1 = origin[1] - (size / 2.0) + (grid_y * leaf_width);
      y2 = origin[1] - (size / 2.0) + ((grid_y+1) * leaf_width);
      z1 = origin[2] - (size / 2.0) + (grid_z * leaf_width);
      z2 = origin[2] - (size / 2.0) + ((grid_z+1) * leaf_width);
      
      // Determine cut function values at probe points
      probe_f[0] = phi( VectorObj(x1, y1, z1) ); 
      probe_f[1] = phi( VectorObj(x2, y1, z1) ); 
      probe_f[2] = phi( VectorObj(x1, y1, z2) ); 
      probe_f[3] = phi( VectorObj(x2, y1, z2) ); 
      probe_f[4] = phi( VectorObj(x1, y2, z1) ); 
      probe_f[5] = phi( VectorObj(x2, y2, z1) ); 
      probe_f[6] = phi( VectorObj(x1, y2, z2) ); 
      probe_f[7] = phi( VectorObj(x2, y2, z2) ); 
      probe_f[8] = phi( VectorObj((x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0) ); 
      
      // Determine if it is on the surface:
      // Does it contain a non-negative and non-positive vertex?  Zero counts as both.
      non_negative_found = false;
      non_positive_found = false;
      probe_iter = probe_f.begin();
      for ( ; probe_iter != probe_f.end(); ++probe_iter )
      {
         if ( *probe_iter >= 0.0 ) {
            non_negative_found = true;
         }
         
         if ( *probe_iter <= 0.0 ) {
            non_positive_found = true;
         }
         
         if ( non_negative_found && non_positive_found ) {
            break;
         }
      }
      
      // If this octant is on the surface
      if ( non_negative_found && non_positive_found )
      {
         // Create this leaf octant in the octree
         m_octree.createOctant( depth, grid_x, grid_y, grid_z, probe_f, true );

         // Add left adjacent octant to pending list
         if ( grid_x > 0 && 
              visited.find( m_octree.computeOctantHashKey( depth, grid_x-1, grid_y, grid_z ) ) == visited.end() )
         {
            pending.push( m_octree.computeOctantHashKey( depth, grid_x-1, grid_y, grid_z ) );
         }

         // Add right adjacent octant to pending list
         if ( grid_x < grid_size-1 && 
              visited.find( m_octree.computeOctantHashKey( depth, grid_x+1, grid_y, grid_z ) ) == visited.end() )
         {
            pending.push( m_octree.computeOctantHashKey( depth, grid_x+1, grid_y, grid_z ) );
         }

         // Add bottom adjacent octant to pending list
         if ( grid_y > 0 && 
              visited.find( m_octree.computeOctantHashKey( depth, grid_x, grid_y-1, grid_z ) ) == visited.end() )
         {
            pending.push( m_octree.computeOctantHashKey( depth, grid_x, grid_y-1, grid_z ) );
         }

         // Add top adjacent octant to pending list
         if ( grid_y < grid_size-1 && 
              visited.find( m_octree.computeOctantHashKey( depth, grid_x, grid_y+1, grid_z ) ) == visited.end() )
         {
            pending.push( m_octree.computeOctantHashKey( depth, grid_x, grid_y+1, grid_z ) );
         }

         // Add back adjacent octant to pending list
         if ( grid_z > 0 && 
              visited.find( m_octree.computeOctantHashKey( depth, grid_x, grid_y, grid_z-1 ) ) == visited.end() )
         {
            pending.push( m_octree.computeOctantHashKey( depth, grid_x, grid_y, grid_z-1 ) );
         }

         // Add front adjacent octant to pending list
         if ( grid_z < grid_size-1 && 
              visited.find( m_octree.computeOctantHashKey( depth, grid_x, grid_y, grid_z+1 ) ) == visited.end() )
         {
            pending.push( m_octree.computeOctantHashKey( depth, grid_x, grid_y, grid_z+1 ) );
         }
      }
   }
         
   return true;
}

//-------------------------------------------------------------------------
void
Mesh::enforceContinuationCondition(
    std::function<double (const VectorObj& coords)> phi
    )

{
   // Enforce the continuation condition - ensure that only BCC tetrahedra
   // interesect the surface.

   // Get the octree parameters
   double depth = m_octree.getDepth();
   double leaf_width = m_octree.getLeafWidth();
   const VectorObj& origin = m_octree.getOrigin();
   double grid_size = pow( 2.0, depth );
   double size = leaf_width * grid_size;
   
   // Create a set of hash keys that are pending creation
   std::set<long> pending_creation;
   
   // Loop over leaf octants
   OctantHashMap& leaf_octants = m_octree.getOctantHashMap( depth );
   OctantHashMapIterator octant_iter = leaf_octants.begin();
   for ( ; octant_iter != leaf_octants.end(); ++octant_iter )
   {
      Octant* octant = octant_iter->second;
      
      // Get octant center vertex
      Vertex* octant_center = octant->getVertex( 8 );
      
      // Get the grid coordinates
      const VectorObj& grid_coords = octant->getGridCoords();

      // Loop over all faces
      for ( int face = 0; face < 6; ++face )
      {
         // If the adjacent octant already exists, skip this face.
         if ( m_octree.getAdjacentOctant( 
                 depth, grid_size, grid_coords[0], grid_coords[1], grid_coords[2], face ) ) {
            continue;
         }
         
         // If on the border or the adjacent octant will be created, skip this face.
         if ( face == 0 && (grid_coords[0] <= 0 || pending_creation.find( 
                 m_octree.computeOctantHashKey(depth, grid_coords[0]-1, grid_coords[1], grid_coords[2]) ) 
                    != pending_creation.end()) ) {
            continue;
         }    
         else if ( face == 1 && (grid_coords[0] >= grid_size-1 || pending_creation.find( 
                      m_octree.computeOctantHashKey(depth, grid_coords[0]+1, grid_coords[1], grid_coords[2]) ) 
                         != pending_creation.end()) ) {
            continue;
         }
         else if ( face == 2 && (grid_coords[2] <= 0 || pending_creation.find( 
                      m_octree.computeOctantHashKey(depth, grid_coords[0], grid_coords[1], grid_coords[2]-1) ) 
                         != pending_creation.end()) ) {
            continue;
         }
         else if ( face == 3 && (grid_coords[2] >= grid_size-1 || pending_creation.find( 
                      m_octree.computeOctantHashKey(depth, grid_coords[0], grid_coords[1], grid_coords[2]+1) ) 
                         != pending_creation.end()) ) {
            continue;
         }
         else if ( face == 4 && (grid_coords[1] <= 0 || pending_creation.find( 
                      m_octree.computeOctantHashKey(depth, grid_coords[0], grid_coords[1]-1, grid_coords[2]) ) 
                         != pending_creation.end()) ) {
            continue;
         }
         else if ( face == 5 && (grid_coords[1] >= grid_size-1 || pending_creation.find( 
                      m_octree.computeOctantHashKey(depth, grid_coords[0], grid_coords[1]+1, grid_coords[2]) ) 
                         != pending_creation.end()) ) {
            continue;
         }
            
         // Get the face
         Vertex* p1 = NULL;
         Vertex* p2 = NULL;
         Vertex* p3 = NULL;
         Vertex* p4 = NULL;
         Vertex* face_center = NULL;
         m_octree.getOctantFace( 
            depth, grid_size, grid_coords[0], grid_coords[1], grid_coords[2],
            face, p1, p2, p3, p4, face_center );
         
         // If a face has at least one non-positive vertex and one
         // non-negative vertex (zero counts as both), create a leaf
         // octant on the adjoining side of the face.
         // Moreover, if an octant has a corner vertex whose sign is opposite the center vertex,
         // or if either is zero, then create the three leaf octants incident to the corner vertex
         // that share a square face with the octant.  This is handled through the ORS below
         // that will check all 4 vertices of the face for this condition.
         if ( ((p1->getF() <= 0.0 || p2->getF() <= 0.0 || p3->getF() <= 0.0 || p4->getF() <= 0.0) &&
               (p1->getF() >= 0.0 || p2->getF() >= 0.0 || p3->getF() >= 0.0 || p4->getF() >= 0.0)) ||
              p1->getF() * octant_center->getF() <= 0.0 ||
              p2->getF() * octant_center->getF() <= 0.0 ||
              p3->getF() * octant_center->getF() <= 0.0 ||
              p4->getF() * octant_center->getF() <= 0.0 ) 
         {
            switch (face)
            {
               case 0:
                  pending_creation.insert( 
                     m_octree.computeOctantHashKey(depth, grid_coords[0]-1, grid_coords[1], grid_coords[2]) );
                  break;
               case 1:
                  pending_creation.insert( 
                     m_octree.computeOctantHashKey(depth, grid_coords[0]+1, grid_coords[1], grid_coords[2]) );
                  break;
               case 2:
                  pending_creation.insert( 
                     m_octree.computeOctantHashKey(depth, grid_coords[0], grid_coords[1], grid_coords[2]-1) );
                  break;
               case 3:
                  pending_creation.insert( 
                     m_octree.computeOctantHashKey(depth, grid_coords[0], grid_coords[1], grid_coords[2]+1) );
                  break;
               case 4:
                  pending_creation.insert( 
                     m_octree.computeOctantHashKey(depth, grid_coords[0], grid_coords[1]-1, grid_coords[2]) );
                  break;
               case 5:
                  pending_creation.insert( 
                     m_octree.computeOctantHashKey(depth, grid_coords[0], grid_coords[1]+1, grid_coords[2]) );
                  break;
            };
         } 
      }    
   }     
   
   // Create the octants that were pending creation
   int grid_x = 0;
   int grid_y = 0;
   int grid_z = 0;
   std::vector<double> probe_f;
   probe_f.resize( 9 );
   std::set<long>::iterator create_iter = pending_creation.begin();
   for ( ; create_iter != pending_creation.end(); ++create_iter )
   {
      // Get the grid coordinates from the hash key
      m_octree.inverseOctantHashKey( depth, *create_iter, grid_x, grid_y, grid_z );
      
      // Determine 3D coordinate range of probe points for current octant
      double x1 = origin[0] - (size / 2.0) + (grid_x * leaf_width);
      double x2 = origin[0] - (size / 2.0) + ((grid_x+1) * leaf_width);
      double y1 = origin[1] - (size / 2.0) + (grid_y * leaf_width);
      double y2 = origin[1] - (size / 2.0) + ((grid_y+1) * leaf_width);
      double z1 = origin[2] - (size / 2.0) + (grid_z * leaf_width);
      double z2 = origin[2] - (size / 2.0) + ((grid_z+1) * leaf_width);
      
      // Determine cut function values at probe points
      probe_f[0] = phi( VectorObj(x1, y1, z1) ); 
      probe_f[1] = phi( VectorObj(x2, y1, z1) ); 
      probe_f[2] = phi( VectorObj(x1, y1, z2) ); 
      probe_f[3] = phi( VectorObj(x2, y1, z2) ); 
      probe_f[4] = phi( VectorObj(x1, y2, z1) ); 
      probe_f[5] = phi( VectorObj(x2, y2, z1) ); 
      probe_f[6] = phi( VectorObj(x1, y2, z2) ); 
      probe_f[7] = phi( VectorObj(x2, y2, z2) ); 
      probe_f[8] = phi( VectorObj((x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0) );
      
      // Create the leaf octant
      m_octree.createOctant( depth, grid_x, grid_y, grid_z, probe_f, true );      
   }
}

//-----------------------------------------------------
void
Mesh::buildOctree()
{
   // Build the octree - create ancestors of leaf octants
   
   // Get octree parameters
   int depth = m_octree.getDepth();
   double leaf_width = m_octree.getLeafWidth();
      
   // Iterate over depths leading to the root (0 depth)
   Octant* child0 = NULL;
   Octant* child1 = NULL;
   Octant* child2 = NULL;
   Octant* child3 = NULL;
   Octant* child4 = NULL;
   Octant* child5 = NULL;
   Octant* child6 = NULL;
   Octant* child7 = NULL;
   std::vector<double> probe_f;
   for ( int current_depth = depth-1; current_depth >= 0; --current_depth )
   {
      // Get current grid size, spatial size, and octant width
      int grid_size = (int) pow( 2.0, (double) current_depth );
   
      // Loop over all posible octants
      for ( int grid_x = 0; grid_x < grid_size; ++grid_x )
      {
         for ( int grid_y = 0; grid_y < grid_size; ++grid_y )
         {
            for ( int grid_z = 0; grid_z < grid_size; ++grid_z )
            {
               // Attempt to retrieve 8 children
               child0 = m_octree.getOctant( current_depth+1, grid_x*2,   grid_y*2,   grid_z*2 );
               child1 = m_octree.getOctant( current_depth+1, grid_x*2+1, grid_y*2,   grid_z*2 );
               child2 = m_octree.getOctant( current_depth+1, grid_x*2,   grid_y*2,   grid_z*2+1 );
               child3 = m_octree.getOctant( current_depth+1, grid_x*2+1, grid_y*2,   grid_z*2+1 );
               child4 = m_octree.getOctant( current_depth+1, grid_x*2,   grid_y*2+1, grid_z*2 );
               child5 = m_octree.getOctant( current_depth+1, grid_x*2+1, grid_y*2+1, grid_z*2 );
               child6 = m_octree.getOctant( current_depth+1, grid_x*2,   grid_y*2+1, grid_z*2+1 );
               child7 = m_octree.getOctant( current_depth+1, grid_x*2+1, grid_y*2+1, grid_z*2+1 );

               // If there is at least one child, create this octant
               if ( child0 != NULL || child1 != NULL || child2 != NULL || child3 != NULL ||
                    child4 != NULL || child5 != NULL || child6 != NULL || child7 != NULL )
               {
                  // (Note: we will hold off on calculating cut function values for the vertices
                  // until we need them during the background grid generation).
                  m_octree.createOctant( current_depth, grid_x, grid_y, grid_z, probe_f, false );
               }
            }
         }
      }
   }
}

//-----------------------------------------------------
void
Mesh::enforceWeakBalanceCondition()
{
   // Enforce the weak balance condition - make sure octants span by no more than a factor
   // of 2.
   
   // Get octree parameters
   int depth = m_octree.getDepth();
   
   // If the tree has depth less or equal to 2 levels, no need for this.
   if ( depth <= 1 ) {
      return;
   }
   
   // Loop over all depths 1-away from the leaves, to the root.
   Octant* octant = NULL;
   std::vector<double> probe_f;
   for ( int current_depth = depth - 1; current_depth >= 0; --current_depth )
   {
      // Get current grid size, spatial size, and octant width
      int grid_size = (int) pow( 2.0, (double) current_depth );
      
      // Loop over all posible octants
      for ( int grid_x = 0; grid_x < grid_size; ++grid_x )
      {
         for ( int grid_y = 0; grid_y < grid_size; ++grid_y )
         {
            for ( int grid_z = 0; grid_z < grid_size; ++grid_z )
            {
               // Skip if the octant already exists
               octant = m_octree.getOctant( current_depth, grid_x, grid_y, grid_z );
               if ( octant != NULL ) {
                  continue;
               }
               
               // Determine the vertex grid coordinates for this potential octant
               GridCoordVector vertex_grid_coords;
               vertex_grid_coords.resize( 9 );
               m_octree.computeVertexGridCoordinates( 
                  current_depth, grid_x, grid_y, grid_z, vertex_grid_coords ); 
            
               // If any of the 12 edges of this octant contains a mid point, 
               // it means that an edge half of its length intersects this octant.              
               if ( m_octree.getVertex( (vertex_grid_coords[0] + vertex_grid_coords[1]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[0] + vertex_grid_coords[2]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[1] + vertex_grid_coords[3]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[2] + vertex_grid_coords[3]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[4] + vertex_grid_coords[5]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[4] + vertex_grid_coords[6]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[5] + vertex_grid_coords[7]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[6] + vertex_grid_coords[7]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[2] + vertex_grid_coords[6]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[3] + vertex_grid_coords[7]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[0] + vertex_grid_coords[4]) / 2.0 ) != NULL ||
                    m_octree.getVertex( (vertex_grid_coords[1] + vertex_grid_coords[5]) / 2.0 ) != NULL )
               {
                  // Create this octant and all its ancestors                     
                  for ( int d = current_depth; d >= 0; --d )
                  {
                     int divide_factor = (int) pow( 2.0, current_depth-d );
                  
                     if ( m_octree.getOctant( 
                        d, grid_x / divide_factor, grid_y / divide_factor, grid_z / divide_factor ) == NULL )
                     {
                        // Note: We are holding off on computing cut function values until we need them later.
                        m_octree.createOctant( 
                           d, grid_x / divide_factor, grid_y / divide_factor, grid_z / divide_factor, probe_f, false );
                     }
                  }
               }
            }
         }
      }      
   }
}

//-----------------------------------------------------
void
Mesh::computeBackgroundGrid( 
    std::function<double (const VectorObj& coords)> phi,
   VertexList& vertex_subset, TetrahedronVector& background_grid )
{
   // Convert the balanced octree to a background grid of tetrahedra to be stuffed
   // later.  This will return a subset of vertices in the background grid,
   // and the background grid tetrahedra.

   // Clear return lists
   vertex_subset.clear();
   background_grid.clear();   
   
   // Get octree parameters
   int depth = m_octree.getDepth();
   
   // Loop over all depths, from the leaves to the root
   for ( int current_depth = depth; current_depth >= 0; --current_depth )
   {    
      // Create a set of visited octants (per depth) to avoid duplication of 
      // tetrahedra that span two octants.
      std::set<long> visited;
   
      // Calculate grid size at current depth
      int grid_size = (int) pow( 2.0, (double) current_depth );
   
      // Loop over octants at current depth
      OctantHashMap& octants = m_octree.getOctantHashMap( current_depth );
      OctantHashMapIterator octant_iter = octants.begin();
      for ( ; octant_iter != octants.end(); ++octant_iter )
      {      
         // o = current octant
         Octant* o = octant_iter->second; 

         // c = center vertex of octant
         Vertex* c = o->getVertex(8);

         // Skip if this octant is not a leaf and has a negative center 
         // point (i.e. octant is outside of surface)
         if ( current_depth != depth )
         { 
            // Make sure to compute the cut function value, if not already.
            if ( c->getFComputed() == false ) {
               c->setF( phi( c->getCoords() ) );   
            }
         
            if ( c->getF() < 0.0 ) {
               continue;
            }
         }
                  
         // Get the grid coordinates of the octant
         const VectorObj& grid_coords = o->getGridCoords();
         
         // For each square face s of octant
         for ( int s = 0; s < 6; ++s )
         {         
            // Get the face (corners and center vertex d)
            Vertex* face_corner1 = NULL; 
            Vertex* face_corner2 = NULL; 
            Vertex* face_corner3 = NULL; 
            Vertex* face_corner4 = NULL; 
            Vertex* d = NULL; 
            m_octree.getOctantFace( 
               current_depth, grid_size, grid_coords[0], grid_coords[1], grid_coords[2],
               s, face_corner1, face_corner2, face_corner3, face_corner4, d );
           
            // If there is no vertex at the center of s
            if ( d == NULL )
            {
               // If s is shared with another octant o_prime the same size as o
               Octant* o_prime = m_octree.getAdjacentOctant( 
                  current_depth, grid_size, grid_coords[0], grid_coords[1], grid_coords[2], s );
               if ( o_prime != NULL )                     
               {
                  // If we haven't visited o_prime yet (avoid duplication)
                  if ( visited.find( m_octree.computeOctantHashKey(
                          current_depth, o_prime->getGridCoords()) ) == visited.end() )
                  {           
                     // c_prime = the center vertex of o_prime
                     Vertex* c_prime = o_prime->getVertex(8);
                     
                     // For each edge e of the square s
                     EdgeIndexList edge_indices;
                     m_octree.getEdgeIndicesForFace(s, edge_indices);
                     EdgeIndexListIterator edge_iter = edge_indices.begin();
                     for ( ; edge_iter != edge_indices.end(); ++edge_iter )
                     {
                        // e = current edge index
                        int e = *edge_iter;
                     
                        // Get edge end points p1,p2 and center point m
                        Vertex* p1 = NULL; 
                        Vertex* p2 = NULL;
                        Vertex* m = NULL; 
                        m_octree.getEdge( 
                           current_depth, grid_size, grid_coords[0], grid_coords[1], grid_coords[2],
                           e, p1, p2, m );
                           
                        // If there is no vertex at the midpoint of e
                        if ( m == NULL )
                        {                     
                           // Create the BCC tetrahedron conv(p1,p2,c,c_prime)
                           createBackgroundGridTetrahedron(
                              phi, p1, p2, c, c_prime,
                              BLACK_EDGE, RED_EDGE, RED_EDGE, RED_EDGE, RED_EDGE, BLACK_EDGE,
                              BCC_TETRAHEDRON, vertex_subset, background_grid );
                        }
                        else
                        {
                            // Create the bisected BCC tet conv(p1, m, c, c_prime)
                           createBackgroundGridTetrahedron(
                              phi, p1, m, c, c_prime,
                              BLACK_EDGE, RED_EDGE, RED_EDGE, BLUE_EDGE, BLUE_EDGE, BLACK_EDGE,
                              BISECTED_BCC_TETRAHEDRON, vertex_subset, background_grid );
                           
                           // Create the bisected BCC tet conv(p2, m, c, c_prime)
                           createBackgroundGridTetrahedron(
                              phi, p2, m, c, c_prime,
                              BLACK_EDGE, RED_EDGE, RED_EDGE, BLUE_EDGE, BLUE_EDGE, BLACK_EDGE,
                              BISECTED_BCC_TETRAHEDRON, vertex_subset, background_grid );
                        } 
                     }
                  }
               }
               else
               {
                  // s is shared with a larger octant or the boundary
                                    
                  // Create two half pyramids filling the pyramid conv(s union c)
                  // Choose diagonal adjoining corner or center of c's parent.
                  if ( o->getParent()->containsVertex(face_corner1) ||
                       o->getParent()->containsVertex(face_corner3) )
                  {
                     createBackgroundGridTetrahedron(
                        phi, face_corner1, face_corner2, face_corner3, c,
                        BLACK_EDGE, BLUE_EDGE, RED_EDGE, BLACK_EDGE, RED_EDGE, RED_EDGE,
                        HALF_PYRAMID_TETRAHEDRON, vertex_subset, background_grid );                     
                     createBackgroundGridTetrahedron(
                        phi, face_corner1, face_corner4, face_corner3, c,
                        BLACK_EDGE, BLUE_EDGE, RED_EDGE, BLACK_EDGE, RED_EDGE, RED_EDGE,
                        HALF_PYRAMID_TETRAHEDRON, vertex_subset, background_grid );                    
                  }
                  else if ( o->getParent()->containsVertex(face_corner2) ||
                            o->getParent()->containsVertex(face_corner4) )
                  {
                     createBackgroundGridTetrahedron(
                        phi, face_corner2, face_corner3, face_corner4, c,
                        BLACK_EDGE, BLUE_EDGE, RED_EDGE, BLACK_EDGE, RED_EDGE, RED_EDGE,
                        HALF_PYRAMID_TETRAHEDRON, vertex_subset, background_grid );                     
                     createBackgroundGridTetrahedron(
                        phi, face_corner2, face_corner1, face_corner4, c,
                        BLACK_EDGE, BLUE_EDGE, RED_EDGE, BLACK_EDGE, RED_EDGE, RED_EDGE,
                        HALF_PYRAMID_TETRAHEDRON, vertex_subset, background_grid );                     
                  }
               } 
            }
            else
            {
               // There is a vertex d at the center of square s
               
               // For each edge e of the square s
               EdgeIndexList edge_indices;
               m_octree.getEdgeIndicesForFace(s, edge_indices);
               EdgeIndexListIterator edge_iter = edge_indices.begin();
               for ( ; edge_iter != edge_indices.end(); ++edge_iter )
               {
                  // e = current edge index
                  int e = *edge_iter;
                  
                  // Get edge end points p1,p2 and center point m
                  Vertex* p1 = NULL; 
                  Vertex* p2 = NULL;
                  Vertex* m = NULL; 
                  m_octree.getEdge( 
                     current_depth, grid_size, grid_coords[0], grid_coords[1], grid_coords[2],
                     e, p1, p2, m );
                     
                  // If there is no vertex at the midpoint of e
                  if ( m == NULL )
                  {
                     // Create a bisected BCC tetrahedron conv(p1,p2,d,c)
                     createBackgroundGridTetrahedron(
                        phi, p1, p2, d, c,
                        BLACK_EDGE, BLUE_EDGE, RED_EDGE, BLUE_EDGE, RED_EDGE, BLACK_EDGE,
                        BISECTED_BCC_TETRAHEDRON, vertex_subset, background_grid );
                  } 
                  else
                  {
                     // m = midpoint vertex of e
                     
                     // Loop over all children and determine if they contain p1, p2
                     bool child_contains_p1 = false;
                     bool child_contains_p2 = false;
                     OctantVector& children = o->getChildren();
                     OctantVectorIterator child_iter = children.begin();
                     for ( ; child_iter != children.end(); ++child_iter )
                     {
                        if ( *child_iter == NULL ) {
                           continue;
                        }
                     
                        if ( (*child_iter)->containsVertex( p1 ) ) {
                           child_contains_p1 = true;
                        }
                        if ( (*child_iter)->containsVertex( p2 ) ) {
                           child_contains_p2 = true;
                        }
                        
                        if ( child_contains_p1 && child_contains_p2 ) {
                           break;
                        }                        
                     }
                     
                     // If o has no child with vertex p1
                     if ( child_contains_p1 == false )
                     {
                        // Create quadrisected BCC tetrahedron conv(p1,m,d,c)
                        createBackgroundGridTetrahedron(
                           phi, p1, m, d, c,
                           BLACK_EDGE, BLUE_EDGE, RED_EDGE, BLACK_EDGE, BLUE_EDGE, BLACK_EDGE,
                           QUADRISECTED_BCC_TETRAHEDRON, vertex_subset, background_grid );
                     }
                     
                     // If o has no child with vertex p2
                     if ( child_contains_p2 == false )
                     {
                        // Create quadrisected BCC tetrahedron conv(p1,m,d,c)
                        createBackgroundGridTetrahedron(
                           phi, p2, m, d, c,
                           BLACK_EDGE, BLUE_EDGE, RED_EDGE, BLACK_EDGE, BLUE_EDGE, BLACK_EDGE,
                           QUADRISECTED_BCC_TETRAHEDRON, vertex_subset, background_grid );
                     }
                  }                  
               }
            }
         }
         
         // Mark that this octant has been visited, to avoid duplicate tetrahedra
         visited.insert( m_octree.computeOctantHashKey( current_depth, grid_coords ) );
      }
   }
}

//-----------------------------------------------------
void 
Mesh::createBackgroundGridTetrahedron(
    std::function<double (const VectorObj& coords)> phi,
   Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4,
   edge_color_type edge_color1, 
   edge_color_type edge_color2, 
   edge_color_type edge_color3, 
   edge_color_type edge_color4, 
   edge_color_type edge_color5, 
   edge_color_type edge_color6,
   tetrahedron_type type, 
   VertexList& vertex_subset, TetrahedronVector& background_grid )
{
   // Create a background grid tetrahedron with the given 4 vertices.
   // Pass the 6 edge colors for new edges:
   // 1 = Between vertices 1 and 2
   // 2 = Between vertices 1 and 3
   // 3 = Between vertices 1 and 4
   // 4 = Between vertices 2 and 3
   // 5 = Between vertices 2 and 4
   // 6 = Between vertices 3 and 4
   // Pass the tetrahedron type: BCC, half BCC, quad BCC, half pyramid
   // This will add the new tetrahedron to the background grid list,
   // create the edges between the vertices, and update the 
   // vertex subset list appropriately.
   
   // Create new tetrahedron
   Tetrahedron* tetrahedron = new Tetrahedron();
   
   // Set its vertices
   tetrahedron->setVertices( v1, v2, v3, v4 );

   // Set the type
   tetrahedron->setType( type );

   // If any of the vertices were just added to the background grid 
   // for the first time, add it to the vertex subset.
   // Also compute the vertex's cut function value, if it hasn't already.
   if ( v1->getEdgeList().empty() ) 
   {
      v1->setInSubset( true );
      vertex_subset.push_back( v1 );
      
      if ( v1->getFComputed() == false ) {
         v1->setF( phi( v1->getCoords() ) );
      }
   }
   if ( v2->getEdgeList().empty() ) 
   {
      v2->setInSubset( true );
      vertex_subset.push_back( v2 );

      if ( v2->getFComputed() == false ) {
         v2->setF( phi( v2->getCoords() ) );
      }
   }
   if ( v3->getEdgeList().empty() ) 
   {
      v3->setInSubset( true );
      vertex_subset.push_back( v3 );

      if ( v3->getFComputed() == false ) {
         v3->setF( phi( v3->getCoords() ) );
      }
   }
   if ( v4->getEdgeList().empty() ) 
   {
      v4->setInSubset( true );
      vertex_subset.push_back( v4 );

      if ( v4->getFComputed() == false ) {
         v4->setF( phi( v4->getCoords() ) );
      }
   }
         
   // Create the 6 edges
   m_octree.createEdge( v1, v2, edge_color1 );
   m_octree.createEdge( v1, v3, edge_color2 );
   m_octree.createEdge( v1, v4, edge_color3 );
   m_octree.createEdge( v2, v3, edge_color4 );
   m_octree.createEdge( v2, v4, edge_color5 );
   m_octree.createEdge( v3, v4, edge_color6 );   
   
   // Add it to the background grid list
   background_grid.push_back( tetrahedron ); 
}
