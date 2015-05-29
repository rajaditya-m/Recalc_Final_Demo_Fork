
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

#ifndef _TETRAHEDRON_H_
#define _TETRAHEDRON_H_

#include "Common.h"
#include "Edge.h"
#include "Util/vectorObj.h"

// Forward Class Declarations
class Vertex;

// Tetrahedron class
class Tetrahedron
{
public:
  
   Tetrahedron();
   virtual ~Tetrahedron();

   void
   getVertices(
      Vertex*& vertex1,
      Vertex*& vertex2,
      Vertex*& vertex3,
      Vertex*& vertex4 ) const;
   void
   setVertices(
      Vertex* vertex1,
      Vertex* vertex2,
      Vertex* vertex3,
      Vertex* vertex4 );
   
   tetrahedron_type
   getType() const;
   void
   setType( tetrahedron_type type );
      
   int
   getPositiveCount() const;
   int
   getNegativeCount() const;
   int
   getZeroCount() const;

   Edge*
   getEdge(
      Vertex* start_vertex,
      edge_color_type edge_color,
      double end_f,
      Vertex*& end_vertex ) const;
   
   void   
   connectVertices( VertexVector& mesh_vertices );
   
   void
   getOtherVertices(
      Vertex* vertex, 
      Vertex*& vertex_a,
      Vertex*& vertex_b,
      Vertex*& vertex_c ) const;

   Vertex*
   getOtherVertex(
      Vertex* vertex_a,
      Vertex* vertex_b,
      Vertex* vertex_c ) const;
   
   Vertex*
   getClosestVertex( const VectorObj& position ) const;
   
   bool
   isWellConditioned( 
      double min_dihedral_angle,
      double tetrahedron_tolerance ) const;
   
   bool
   containsVertex( Vertex* vertex ) const;
         
private:

   Vertex* m_vertex1;            // Four vertices
   Vertex* m_vertex2;
   Vertex* m_vertex3;
   Vertex* m_vertex4;

   tetrahedron_type m_type;      // BCC, Bisected BCC, Quadrisected BCC, Half Pyramid
   
};

#endif 
