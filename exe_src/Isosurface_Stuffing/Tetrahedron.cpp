
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

#include "Tetrahedron.h"
#include "Vertex.h"

//-----------------------------------------------------
Tetrahedron::Tetrahedron()
:  m_vertex1 (NULL),
   m_vertex2 (NULL),
   m_vertex3 (NULL),
   m_vertex4 (NULL),
   m_type (BCC_TETRAHEDRON)
{
   // Constructor
}

//-----------------------------------------------------
Tetrahedron::~Tetrahedron()
{
   // Destructor
}

//-----------------------------------------------------
void
Tetrahedron::getVertices(
   Vertex*& vertex1,
   Vertex*& vertex2,
   Vertex*& vertex3,
   Vertex*& vertex4 ) const
{
   vertex1 = m_vertex1;
   vertex2 = m_vertex2;
   vertex3 = m_vertex3;
   vertex4 = m_vertex4;
}

//-----------------------------------------------------
void
Tetrahedron::setVertices(
   Vertex* vertex1,
   Vertex* vertex2,
   Vertex* vertex3,
   Vertex* vertex4 )
{
   m_vertex1 = vertex1;
   m_vertex2 = vertex2;
   m_vertex3 = vertex3;
   m_vertex4 = vertex4;
}

//-----------------------------------------------------
tetrahedron_type
Tetrahedron::getType() const
{
   return m_type;
}
   
//-----------------------------------------------------
void
Tetrahedron::setType( tetrahedron_type type )
{
   m_type = type;
}

//-----------------------------------------------------
int
Tetrahedron::getPositiveCount() const
{
   int count = 0;
   
   if ( m_vertex1 != NULL && m_vertex1->getF() > 0.0 ) {
      count++;
   }
   if ( m_vertex2 != NULL && m_vertex2->getF() > 0.0 ) {
      count++;
   }
   if ( m_vertex3 != NULL && m_vertex3->getF() > 0.0 ) {
      count++;
   }
   if ( m_vertex4 != NULL && m_vertex4->getF() > 0.0 ) {
      count++;
   }
   
   return count;
}

//-----------------------------------------------------
int
Tetrahedron::getNegativeCount() const
{
   int count = 0;
   
   if ( m_vertex1 != NULL && m_vertex1->getF() < 0.0 ) {
      count++;
   }
   if ( m_vertex2 != NULL && m_vertex2->getF() < 0.0 ) {
      count++;
   }
   if ( m_vertex3 != NULL && m_vertex3->getF() < 0.0 ) {
      count++;
   }
   if ( m_vertex4 != NULL && m_vertex4->getF() < 0.0 ) {
      count++;
   }
   
   return count;
}

//-----------------------------------------------------
int
Tetrahedron::getZeroCount() const
{
   int count = 0;
   
   if ( m_vertex1 != NULL && m_vertex1->getF() == 0.0 ) {
      count++;
   }
   if ( m_vertex2 != NULL && m_vertex2->getF() == 0.0 ) {
      count++;
   }
   if ( m_vertex3 != NULL && m_vertex3->getF() == 0.0 ) {
      count++;
   }
   if ( m_vertex4 != NULL && m_vertex4->getF() == 0.0 ) {
      count++;
   }
   
   return count;
}

//-----------------------------------------------------
Edge*
Tetrahedron::getEdge(
   Vertex* start_vertex,
   edge_color_type edge_color,
   double end_f,
   Vertex*& end_vertex ) const
{
   // Retrieve an edge from a start vertex, an edge color,
   // and a positive/zero/negative end vertex value.
   // Also returns the end vertex.
   
   // Loop through all adjacent edges to the start vertex
   end_vertex = NULL;
   Edge* edge = NULL;
   EdgeList& edge_list = start_vertex->getEdgeList();
   EdgeListIterator edge_iter = edge_list.begin();
   for ( ; edge_iter != edge_list.end(); ++edge_iter )
   {
      edge = *edge_iter;
      
      // Determine end vertex
      if ( edge->getVertex1() == start_vertex ) {
         end_vertex = edge->getVertex2();
      }
      else {
         end_vertex = edge->getVertex1();
      }
      
      // Skip if this edge is not in the tetrahedron
      if ( end_vertex != m_vertex1 && 
           end_vertex != m_vertex2 && 
           end_vertex != m_vertex3 && 
           end_vertex != m_vertex4 )
      {
         continue;
      }
      
      // If color and sign match, we found our edge
      if ( edge->getColor() == edge_color &&
           ( ( end_f < 0.0 && end_vertex->getF() < 0.0 ||
               end_f == 0.0 && end_vertex->getF() == 0.0 ||
               end_f > 0.0 && end_vertex->getF() > 0.0 ) ) )
      {
         return edge;
      }
   }
   
   // Not found
   end_vertex = NULL;
   return NULL;
}

//-----------------------------------------------------
void   
Tetrahedron::connectVertices( VertexVector& mesh_vertices )
{
   // Hook up vertex connections for this tetrahedron
   
   // If this is a newly encountered vertex, add it to the list of
   // mesh vertices
   if ( m_vertex1->getTetrahedronVector().empty() ) {
      mesh_vertices.push_back( m_vertex1 );
   }
   if ( m_vertex2->getTetrahedronVector().empty() ) {
      mesh_vertices.push_back( m_vertex2 );
   }
   if ( m_vertex3->getTetrahedronVector().empty() ) {
      mesh_vertices.push_back( m_vertex3 );
   }
   if ( m_vertex4->getTetrahedronVector().empty() ) {
      mesh_vertices.push_back( m_vertex4 );
   }
            
   // Add this tetrahedron to all the vertices' adjacent tetrahedron
   m_vertex1->getTetrahedronVector().push_back( this );
   m_vertex2->getTetrahedronVector().push_back( this );
   m_vertex3->getTetrahedronVector().push_back( this );
   m_vertex4->getTetrahedronVector().push_back( this );
}

//-----------------------------------------------------
void
Tetrahedron::getOtherVertices(
   Vertex* vertex, 
   Vertex*& vertex_a,
   Vertex*& vertex_b,
   Vertex*& vertex_c ) const
{
   // Return the other three vertices to the given one.
   
   if ( m_vertex1 == vertex ) {
      vertex_a = m_vertex2;
      vertex_b = m_vertex3;
      vertex_c = m_vertex4;
   }
   else if ( m_vertex2 == vertex ) {
      vertex_a = m_vertex1;
      vertex_b = m_vertex3;
      vertex_c = m_vertex4;
   }
   else if ( m_vertex3 == vertex ) {
      vertex_a = m_vertex1;
      vertex_b = m_vertex2;
      vertex_c = m_vertex4;
   }
   else {
      vertex_a = m_vertex1;
      vertex_b = m_vertex2;
      vertex_c = m_vertex3;
   }
}

//-----------------------------------------------------
Vertex*
Tetrahedron::getOtherVertex(
   Vertex* vertex_a,
   Vertex* vertex_b,
   Vertex* vertex_c ) const
{
   // Get the other vertex to the 3 that are given
   
   if ( (vertex_a == m_vertex1 || vertex_a == m_vertex2 || vertex_a == m_vertex3) &&
        (vertex_b == m_vertex1 || vertex_b == m_vertex2 || vertex_b == m_vertex3) &&
        (vertex_c == m_vertex1 || vertex_c == m_vertex2 || vertex_c == m_vertex3) )
   {
      return m_vertex4;
   }
   else if ( (vertex_a == m_vertex1 || vertex_a == m_vertex2 || vertex_a == m_vertex4) &&
             (vertex_b == m_vertex1 || vertex_b == m_vertex2 || vertex_b == m_vertex4) &&
             (vertex_c == m_vertex1 || vertex_c == m_vertex2 || vertex_c == m_vertex4) )
   {
      return m_vertex3;
   }
   else if ( (vertex_a == m_vertex1 || vertex_a == m_vertex3 || vertex_a == m_vertex4) &&
             (vertex_b == m_vertex1 || vertex_b == m_vertex3 || vertex_b == m_vertex4) &&
             (vertex_c == m_vertex1 || vertex_c == m_vertex3 || vertex_c == m_vertex4) )
   {
      return m_vertex2;
   }
   else
   {
      return m_vertex1;
   }
}

//-----------------------------------------------------
Vertex*
Tetrahedron::getClosestVertex( const VectorObj& position ) const
{
   // Find the closest vertex to the given coordinates
   
   // Calculate distance to each vertex
   double distance1 = (m_vertex1->getCoords() - position).length();
   double distance2 = (m_vertex2->getCoords() - position).length();
   double distance3 = (m_vertex3->getCoords() - position).length();
   double distance4 = (m_vertex4->getCoords() - position).length();
   
   // Return vertex with shortest distance
   if ( distance1 <= distance2 && distance1 <= distance3 && distance1 <= distance4 ) {
      return m_vertex1;
   }   
   if ( distance2 <= distance1 && distance2 <= distance3 && distance2 <= distance4 ) {
      return m_vertex2;
   }   
   if ( distance3 <= distance1 && distance3 <= distance2 && distance3 <= distance4 ) {
      return m_vertex3;
   }
   else {
      return m_vertex4;
   }   
}

//-----------------------------------------------------
bool
Tetrahedron::isWellConditioned( 
   double min_dihedral_angle,
   double tetrahedron_tolerance ) const
{
   // Determine if this tetrahedron has all dihedral angles
   // above the given minimum angle (in radians).
   // Also pass a tolerance for checking very flat tetrahedron.
   
   // Get vertex coordinates
   VectorObj& coords1 = m_vertex1->getCoords();
   VectorObj& coords2 = m_vertex2->getCoords();
   VectorObj& coords3 = m_vertex3->getCoords();
   VectorObj& coords4 = m_vertex4->getCoords();
   
   // Compute normals
   // Face1: V1, V2, V3
   double dot_prod = 0.0;
   VectorObj normal1 = (coords2 - coords1).cross(coords3 - coords1);
   normal1.normalize();
   dot_prod = (coords4 - coords1).dot( normal1 );
   if ( dot_prod > -tetrahedron_tolerance && dot_prod < tetrahedron_tolerance ) {
      return false;
   }
   if ( dot_prod > 0.0 ) {
      normal1 *= -1.0;
   }

   // Face2: V1, V2, V4
   VectorObj normal2 = (coords2 - coords1).cross(coords4 - coords1);
   normal2.normalize();
   dot_prod = (coords3 - coords1).dot( normal2 );
   if ( dot_prod > -tetrahedron_tolerance && dot_prod < tetrahedron_tolerance ) {
      return false;
   }
   if ( dot_prod > 0.0 ) {
      normal2 *= -1.0;
   }

   // Face3: V1, V3, V4
   VectorObj normal3 = (coords3 - coords1).cross(coords4 - coords1);
   normal3.normalize();
   dot_prod = (coords2 - coords1).dot( normal3 );
   if ( dot_prod > -tetrahedron_tolerance && dot_prod < tetrahedron_tolerance ) {
      return false;
   }
   if ( dot_prod > 0.0 ) {
      normal3 *= -1.0;
   }

   // Face4: V2, V3, V4
   VectorObj normal4 = (coords3 - coords2).cross(coords4 - coords2);
   normal4.normalize();
   dot_prod = (coords1 - coords2).dot( normal4 );
   if ( dot_prod > -tetrahedron_tolerance && dot_prod < tetrahedron_tolerance ) {
      return false;
   }
   if ( dot_prod > 0.0 ) {
      normal4 *= -1.0;
   }
   
   // Check face1 - face2
   double angle = 0.0;
   angle = acos( (-normal1).dot( normal2 ) );
   if ( angle < min_dihedral_angle ) {
      return false;
   }
 
   // Check face1 - face3
   angle = acos( (-normal1).dot( normal3 ) );
   if ( angle < min_dihedral_angle ) {
      return false;
   }

   // Check face1 - face4
   angle = acos( (-normal1).dot( normal4 ) );
   if ( angle < min_dihedral_angle ) {
      return false;
   }

   // Check face2 - face3
   angle = acos( (-normal2).dot( normal3 ) );
   if ( angle < min_dihedral_angle ) {
      return false;
   }

   // Check face2 - face4
   angle = acos( (-normal2).dot( normal4 ) );
   if ( angle < min_dihedral_angle ) {
      return false;
   }

   // Check face3 - face4
   angle = acos( (-normal3).dot( normal4 ) );
   if ( angle < min_dihedral_angle ) {
      return false;
   }
                       
   return true;
}

//-----------------------------------------------------
bool
Tetrahedron::containsVertex( Vertex* vertex ) const
{
   // Returns if this tetrahedron contains the given vertex
   
   return
      vertex == m_vertex1 || 
      vertex == m_vertex2 ||
      vertex == m_vertex3 ||
      vertex == m_vertex4;
}


