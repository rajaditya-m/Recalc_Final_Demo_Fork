
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

#include "Octree.h"
#include "Vertex.h"
#include "Edge.h"

//-----------------------------------------------------
Octree::Octree()
:  m_depth (0),
   m_leaf_width (0.0)
{
   // Constructor
}

//-----------------------------------------------------
Octree::~Octree()
{
   // Destructor
   
   // Destroy the octants
   if ( m_depth > 0 )
   {
      for ( int i = m_depth; i >= 0; --i ) {
         OctantHashMap& hash_map = m_hash_map_vector[i];
         OctantHashMapIterator octant_iter = hash_map.begin();
         for ( ; octant_iter != hash_map.end(); ++octant_iter )
         {
            delete octant_iter->second;
         }
      }
      m_hash_map_vector.clear();
   }
   
   // Destroy the vertices
   VertexHashMapIterator vertex_iter = m_vertices.begin();
   for ( ; vertex_iter != m_vertices.end(); ++vertex_iter )
   {
      delete vertex_iter->second;
   }
   m_vertices.clear();
   
   // Destroy the edges
   EdgeVectorIterator edge_iter = m_edges.begin();
   for ( ; edge_iter != m_edges.end(); ++edge_iter )
   {
      delete *edge_iter;
   }
   m_edges.clear();   
}

//-----------------------------------------------------
int
Octree::getDepth() const
{
   return m_depth;
}

//-----------------------------------------------------
void 
Octree::setDepth( int depth )
{
   m_depth = depth;
   
   // Resize our hash map vector (plus one level for root)
   m_hash_map_vector.resize( depth+1 );
}

//-----------------------------------------------------
double
Octree::getLeafWidth() const
{
   return m_leaf_width;
}

//-----------------------------------------------------
void
Octree::setLeafWidth( double leaf_width )
{
   m_leaf_width = leaf_width;
}

//-----------------------------------------------------
const VectorObj&
Octree::getOrigin() const
{
   return m_origin;
}

//-----------------------------------------------------
void
Octree::setOrigin( const VectorObj& origin )
{
   m_origin = origin;
}

//-----------------------------------------------------
VertexHashMap&
Octree::getVertices()
{
   return m_vertices;
}

//-----------------------------------------------------
EdgeVector&
Octree::getEdges()
{
   return m_edges;
}

//-----------------------------------------------------
Octant*
Octree::getOctant( int depth, int grid_x, int grid_y, int grid_z )
{
   // Get an octant by the given depth and grid coordinates
   // i.e. (0, 0, 0, 0) = root, (1, 0, 0, 0) = child0, (1, 1, 0, 0) = child1, etc.
   
   // Make sure depth is valid
   if ( depth < 0 || depth >= m_hash_map_vector.size() ) {
      return NULL;
   }
   
   // Attempt to find the octant
   OctantHashMapIterator iter = 
      m_hash_map_vector[depth].find( computeOctantHashKey(depth, grid_x, grid_y, grid_z) );

   // If found, return it
   if ( iter != m_hash_map_vector[depth].end() ) {
      return iter->second;   
   }
   else {
      return NULL;
   }
}

//-----------------------------------------------------
Octant*
Octree::getOctant( int depth, long hash_key )
{
   // Get an octant by the given depth and hash key
   
   // Make sure depth is valid
   if ( depth < 0 || depth >= m_hash_map_vector.size() ) {
      return NULL;
   }
   
   // Attempt to find the octant
   OctantHashMapIterator iter = m_hash_map_vector[depth].find( hash_key );

   // If found, return it
   if ( iter != m_hash_map_vector[depth].end() ) {
      return iter->second;   
   }
   else {
      return NULL;
   }
}

//-----------------------------------------------------
OctantHashMapVector&
Octree::getOctantHashMapVector()
{
   return m_hash_map_vector;
}

//-----------------------------------------------------
OctantHashMap&
Octree::getOctantHashMap( int depth )
{
   // Return the octant hash map at the given depth
   
   return m_hash_map_vector[depth];
}

//-----------------------------------------------------
Octant*
Octree::createOctant(
   int depth, int grid_x, int grid_y, int grid_z, 
   std::vector<double>& probe_f, bool probe_f_provided )
{
   // Create an octant at the given depth and grid coordinates.
   // Optionally pass the cut function values for the 9 probe points (vertices)
   // to be copied over.
   // This will properly hook up children to the octant, and if a parent exists,
   // it will hook up to the parent as well. 
   // Any existent vertices in the octree will be copied over.
   // Returns the new octant, or NULL if an error occurred.

   // For the given depth, compute grid size, spatial size, and octant width
   int grid_size = (int) pow( 2.0, (double) depth );
   double size = m_leaf_width * pow( 2.0, m_depth );
   double octant_width = m_leaf_width * pow( 2.0, m_depth - depth );

   // Make sure depth and grid coordinates are valid
   if ( depth < 0 || depth > m_depth ||
        grid_x < 0 || grid_x >= grid_size || 
        grid_y < 0 || grid_y >= grid_size ||
        grid_z < 0 || grid_z >= grid_size )
   {
      return NULL;
   }
   
   // If this octant already exists, return
   if ( getOctant( depth, grid_x, grid_y, grid_z ) != NULL ) {
      return NULL;
   }

   // Create the new octant
   Octant* new_octant = new Octant(); 

   // Set the depth
   new_octant->setDepth( depth );

   // Set the grid coordinates
   new_octant->setGridCoords( VectorObj( grid_x, grid_y, grid_z ) );
      
   // If this is not a leaf
   if ( depth < m_depth )
   {
      // Attempt to retrieve 8 children
      Octant* child0 = getOctant( depth+1, grid_x*2,   grid_y*2,   grid_z*2 );
      Octant* child1 = getOctant( depth+1, grid_x*2+1, grid_y*2,   grid_z*2 );
      Octant* child2 = getOctant( depth+1, grid_x*2,   grid_y*2,   grid_z*2+1 );
      Octant* child3 = getOctant( depth+1, grid_x*2+1, grid_y*2,   grid_z*2+1 );
      Octant* child4 = getOctant( depth+1, grid_x*2,   grid_y*2+1, grid_z*2 );
      Octant* child5 = getOctant( depth+1, grid_x*2+1, grid_y*2+1, grid_z*2 );
      Octant* child6 = getOctant( depth+1, grid_x*2,   grid_y*2+1, grid_z*2+1 );
      Octant* child7 = getOctant( depth+1, grid_x*2+1, grid_y*2+1, grid_z*2+1 );
      
      // Set children
      new_octant->setChild( 0, child0 );
      new_octant->setChild( 1, child1 );
      new_octant->setChild( 2, child2 );
      new_octant->setChild( 3, child3 );
      new_octant->setChild( 4, child4 );
      new_octant->setChild( 5, child5 );
      new_octant->setChild( 6, child6 );
      new_octant->setChild( 7, child7 );

      // Retrieve the children in vector form (if applicable)
      OctantVector& children = new_octant->getChildren();

      // Set the parent in the children
      for ( int i = 0; i < 8; ++i ) 
      {
         if ( children[i] != NULL ) {
            children[i]->setParent( new_octant );
         }
      }
   }

   // If this is not the root
   if ( depth > 0 )
   {
      // Check if a parent exists
      Octant* parent = getOctant( depth-1, grid_x/2, grid_y/2, grid_z/2 );  
      if ( parent != NULL )
      {
         // Hook up this octant to the parent
         new_octant->setParent( parent );
         
         // Hook up the parent to this octant
         if ( grid_x % 2 == 0 && grid_y % 2 == 0 && grid_z % 2 == 0 ) {
            parent->setChild( 0, new_octant );
         }
         else if ( grid_x % 2 == 1 && grid_y % 2 == 0 && grid_z % 2 == 0 ) {
            parent->setChild( 1, new_octant );
         }
         else if ( grid_x % 2 == 0 && grid_y % 2 == 0 && grid_z % 2 == 1 ) {
            parent->setChild( 2, new_octant );
         }
         else if ( grid_x % 2 == 1 && grid_y % 2 == 0 && grid_z % 2 == 1 ) {
            parent->setChild( 3, new_octant );
         }
         else if ( grid_x % 2 == 0 && grid_y % 2 == 1 && grid_z % 2 == 0 ) {
            parent->setChild( 4, new_octant );
         }
         else if ( grid_x % 2 == 1 && grid_y % 2 == 1 && grid_z % 2 == 0 ) {
            parent->setChild( 5, new_octant );
         }
         else if ( grid_x % 2 == 0 && grid_y % 2 == 1 && grid_z % 2 == 1 ) {
            parent->setChild( 6, new_octant );
         }
         else {
            parent->setChild( 7, new_octant );
         }
      }
   }
   
   // Compute vertex grid coordinates
   GridCoordVector vertex_grid_coords;
   vertex_grid_coords.resize( 9 );
   computeVertexGridCoordinates( depth, grid_x, grid_y, grid_z, vertex_grid_coords );

   // Determine 3D coordinate range of vertices for new octant
   double x1 = m_origin[0] - (size / 2.0) + (grid_x * octant_width);
   double x2 = m_origin[0] - (size / 2.0) + ((grid_x+1) * octant_width);
   double y1 = m_origin[1] - (size / 2.0) + (grid_y * octant_width);
   double y2 = m_origin[1] - (size / 2.0) + ((grid_y+1) * octant_width);
   double z1 = m_origin[2] - (size / 2.0) + (grid_z * octant_width);
   double z2 = m_origin[2] - (size / 2.0) + ((grid_z+1) * octant_width);

   // Loop over all vertices
   VertexVector& vertices = new_octant->getVertices();
   for ( int i = 0; i < 9; ++i )
   {
      // If this vertex exists in the octree, set it in the new octant and continue
      VectorObj& vertex_grid_coord = vertex_grid_coords[i];      
      Vertex* octree_vertex = 
         getVertex( vertex_grid_coord[0], vertex_grid_coord[1], vertex_grid_coord[2] );
      if ( octree_vertex != NULL )
      {
         new_octant->setVertex( i, octree_vertex );
         continue;
      }   
      
      // Create a new vertex
      Vertex* new_vertex = new Vertex();
      
      // Set the grid coordinates
      new_vertex->setGridCoords( vertex_grid_coord );

      // Set color - note, this coloring is only important for leaf octants which
      // might become BCC tetrahedra.
      if ( i != 8 ) {
         new_vertex->setColor( BLACK_VERTEX );
      }
      else {
         new_vertex->setColor( RED_VERTEX );
      }
      
      // Set 3D coordinates
      switch (i)
      {
         case 0:
            new_vertex->setCoords( VectorObj( x1, y1, z1 ) );
            break;
         case 1:
            new_vertex->setCoords( VectorObj( x2, y1, z1 ) );
            break;
         case 2:
            new_vertex->setCoords( VectorObj( x1, y1, z2 ) );
            break;
         case 3:
            new_vertex->setCoords( VectorObj( x2, y1, z2 ) );
            break;
         case 4:
            new_vertex->setCoords( VectorObj( x1, y2, z1 ) );
            break;
         case 5:
            new_vertex->setCoords( VectorObj( x2, y2, z1 ) );
            break;
         case 6:
            new_vertex->setCoords( VectorObj( x1, y2, z2 ) );
            break;
         case 7:
            new_vertex->setCoords( VectorObj( x2, y2, z2 ) );
            break;
         case 8:
            new_vertex->setCoords( VectorObj( (x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0 ) );
            break;
         default:
            break;
      };
      
      // Set the cut function value
      if ( probe_f_provided ) {
         new_vertex->setF( probe_f[i] );
      }

      // Add the new vertex to this octant
      vertices[i] = new_vertex;
      
      // Add the new vertex to the octree
      m_vertices.insert( VertexPair( computeVertexHashKey(vertex_grid_coord), new_vertex ) );
   }

   // Add the octant to the octree hash_map
   m_hash_map_vector[depth].insert( 
      OctantPair( computeOctantHashKey(depth, grid_x, grid_y, grid_z), new_octant ) );   

   return new_octant;
}

//-----------------------------------------------------
long
Octree::computeOctantHashKey( int depth, int grid_x, int grid_y, int grid_z ) const
{
   // Compute the single hash key value from 3 coordinates
      
   double size = pow(2.0, (double) depth);
      
   return (long) ((grid_z * size * size) + (grid_y * size) + grid_x);
}

//-----------------------------------------------------
long
Octree::computeOctantHashKey( int depth, const VectorObj& grid_coords ) const
{
   // Helper method
   
   return computeOctantHashKey( depth, (int) grid_coords[0], (int) grid_coords[1], (int) grid_coords[2] );
}

//-----------------------------------------------------
void
Octree::inverseOctantHashKey( int depth, long hash_key, int& grid_x, int& grid_y, int& grid_z ) const
{
   // Compute the grid coordinates from the given hash key
   
   int size = (int) pow(2.0, (double) depth);
   
   grid_z = hash_key / (size*size);
   grid_y = (hash_key % (size*size)) / size;
   grid_x = (hash_key % (size*size)) % size;
}

//-----------------------------------------------------
long
Octree::computeVertexHashKey( int grid_x, int grid_y, int grid_z ) const
{
   // Compute the single hash key value from 3 coordinates
      
   double size = pow(2.0, (double) m_depth + 1.0) + 1.0;
      
   return (long) ((grid_z * size * size) + (grid_y * size) + grid_x);
}

//-----------------------------------------------------
long
Octree::computeVertexHashKey( const VectorObj& grid_coords ) const
{
   // Helper method
   
   return computeVertexHashKey( (int) grid_coords[0], (int) grid_coords[1], (int) grid_coords[2] );
}

//-----------------------------------------------------
void
Octree::computeVertexGridCoordinates(
   int depth, int octant_grid_x, int octant_grid_y, int octant_grid_z,
   GridCoordVector& vertex_grid_coords )
{
   // Return the vertex grid coordinates for the 9 vertices of the given
   // octant.  The grid coordinate vector is indexed by the octant's
   // vertex index.  The grid size of the vertices is the lowest depth
   // of the octree + 1 (for leaf center points).

   // Determine scale factor
   double factor = pow( 2.0, (double) m_depth - depth + 1.0 );      

   // Determine vertex 0 coordinates
   VectorObj base_coords = 
      VectorObj( octant_grid_x * factor, octant_grid_y * factor, octant_grid_z * factor );

   vertex_grid_coords[0] = base_coords;
   vertex_grid_coords[1] = 
      VectorObj( base_coords[0]+factor, base_coords[1], base_coords[2] ); 
   vertex_grid_coords[2] = 
      VectorObj( base_coords[0], base_coords[1], base_coords[2]+factor ); 
   vertex_grid_coords[3] = 
      VectorObj( base_coords[0]+factor, base_coords[1], base_coords[2]+factor ); 
   vertex_grid_coords[4] = 
      VectorObj( base_coords[0], base_coords[1]+factor, base_coords[2] ); 
   vertex_grid_coords[5] = 
      VectorObj( base_coords[0]+factor, base_coords[1]+factor, base_coords[2] ); 
   vertex_grid_coords[6] = 
      VectorObj( base_coords[0], base_coords[1]+factor, base_coords[2]+factor ); 
   vertex_grid_coords[7] = 
      VectorObj( base_coords[0]+factor, base_coords[1]+factor, base_coords[2]+factor ); 
   vertex_grid_coords[8] = 
      VectorObj( base_coords[0]+(factor/2.0), base_coords[1]+(factor/2.0), base_coords[2]+(factor/2.0) ); 
}

//-----------------------------------------------------
bool
Octree::existSmallerAdjacentOctants( int depth, int grid_size, int grid_x, int grid_y, int grid_z )
{
   // Helper method to check all faces
   
   return existSmallerAdjacentOctants( depth, grid_size, grid_x, grid_y, grid_z, -1 );
}

//-----------------------------------------------------
bool
Octree::existSmallerAdjacentOctants( int depth, int grid_size, int grid_x, int grid_y, int grid_z, int face )
{
   // Do there exist any adjacent octants to the given octant that are smaller 
   // (by exactly a factor of 2).  This will check for adjacent octants to the 
   // given face [0=left, 1=right, 2=back, 3=front, 4=bottom, 5=top, -1=all faces].
   // Note: the target octant does not have to exist for this query.
   // Pass the grid size for efficiency purposes.
      
   // Check the left child
   if ( (face == -1 || face == 0) && grid_x > 0 ) { 
      Octant* left_octant = getOctant( depth, grid_x-1, grid_y, grid_z );
      if ( left_octant != NULL ) 
      {
         // Return true if there are adjacent children
         if ( left_octant->getChild(1) != NULL || 
              left_octant->getChild(3) != NULL ||
              left_octant->getChild(5) != NULL ||
              left_octant->getChild(7) != NULL ) {
            return true;
         }
      }
   }

   // Check the right child
   if ( (face == -1 || face == 1) && grid_x < grid_size-1 ) { 
      Octant* right_octant = getOctant( depth, grid_x+1, grid_y, grid_z );
      if ( right_octant != NULL ) 
      {
         // Return true if there are adjacent children
         if ( right_octant->getChild(0) != NULL || 
              right_octant->getChild(2) != NULL ||
              right_octant->getChild(4) != NULL ||
              right_octant->getChild(6) != NULL ) {
            return true;
         }
      }
   }   

   // Check the back child
   if ( (face == -1 || face == 2) && grid_z > 0 ) { 
      Octant* back_octant = getOctant( depth, grid_x, grid_y, grid_z-1 );
      if ( back_octant != NULL ) 
      {
         // Return true if there are adjacent children
         if ( back_octant->getChild(2) != NULL || 
              back_octant->getChild(3) != NULL ||
              back_octant->getChild(6) != NULL ||
              back_octant->getChild(7) != NULL ) {
            return true;
         }
      }
   }

   // Check the front child
   if ( (face == -1 || face == 3) && grid_z < grid_size-1 ) { 
      Octant* front_octant = getOctant( depth, grid_x, grid_y, grid_z+1 );
      if ( front_octant != NULL ) 
      {
         // Return true if there are adjacent children
         if ( front_octant->getChild(0) != NULL || 
              front_octant->getChild(1) != NULL ||
              front_octant->getChild(4) != NULL ||
              front_octant->getChild(5) != NULL ) {
            return true;
         }
      }
   }

   // Check the bottom child
   if ( (face == -1 || face == 4) && grid_y > 0 ) { 
      Octant* bottom_octant = getOctant( depth, grid_x, grid_y-1, grid_z );
      if ( bottom_octant != NULL ) 
      {
         // Return true if there are adjacent children
         if ( bottom_octant->getChild(4) != NULL || 
              bottom_octant->getChild(5) != NULL ||
              bottom_octant->getChild(6) != NULL ||
              bottom_octant->getChild(7) != NULL ) {
            return true;
         }
      }
   }

   // Check the top child
   if ( (face == -1 || face == 5) && grid_y < grid_size-1 ) { 
      Octant* top_octant = getOctant( depth, grid_x, grid_y+1, grid_z );
      if ( top_octant != NULL ) 
      {
         // Return true if there are adjacent children
         if ( top_octant->getChild(0) != NULL || 
              top_octant->getChild(1) != NULL ||
              top_octant->getChild(2) != NULL ||
              top_octant->getChild(3) != NULL ) {
            return true;
         }
      }
   }
   
   // No smaller adjacent children found
   return false;
}

//-----------------------------------------------------
Octant*
Octree::getAdjacentOctant( int depth, int grid_size, int grid_x, int grid_y, int grid_z, int face )
{
   // Get the adjacent octant of the same size, to the given octant's face, if it exists.
   // Pass the face value [0=left, 1=right, 2=back, 3=front, 4=bottom, 5=top]
   
   // Left child
   if ( face == 0 && grid_x > 0 ) { 
      return getOctant( depth, grid_x-1, grid_y, grid_z );
   }
   else if ( face == 1 && grid_x < grid_size-1 ) {
      return getOctant( depth, grid_x+1, grid_y, grid_z );
   }
   else if ( face == 2 && grid_z > 0 ) { 
      return getOctant( depth, grid_x, grid_y, grid_z-1 );
   }
   else if ( face == 3 && grid_z < grid_size-1 ) {
      return getOctant( depth, grid_x, grid_y, grid_z+1 );
   }   
   else if ( face == 4 && grid_y > 0 ) { 
      return getOctant( depth, grid_x, grid_y-1, grid_z );
   }
   else if ( face == 5 && grid_y < grid_size-1 ) {
      return getOctant( depth, grid_x, grid_y+1, grid_z );
   }

   // Octant is on the border
   return NULL;
}

//-----------------------------------------------------
Vertex*
Octree::getVertex( int grid_x, int grid_y, int grid_z )
{
   // Get the vertex in the octree with the grid coordinates
   
   // Attempt to find the vertex
   VertexHashMapIterator iter = 
      m_vertices.find( computeVertexHashKey(grid_x, grid_y, grid_z) );

   // If found, return it
   if ( iter != m_vertices.end() ) {
      return iter->second;   
   }
   else {
      return NULL;
   }
}

//-----------------------------------------------------
Vertex*
Octree::getVertex( const VectorObj& grid_coords )
{
   // Helper method
   
   return getVertex( (int) grid_coords[0], (int) grid_coords[1], (int) grid_coords[2] );
}

//-----------------------------------------------------
void
Octree::getEdge( 
   int depth, int grid_size, int grid_x, int grid_y, int grid_z, 
   int edge, Vertex*& p1, Vertex*& p2, Vertex*& center )
{
  (void) grid_size;
   // Get an edge of the given octant.
   // Returns the two end points of the edge, and its center point
   // if it exists (NULL otherwise). 
   // Pass an edge value of:
   // 0 = Between vertices 0 and 1
   // 1 = Between vertices 1 and 3
   // 2 = Between vertices 0 and 2
   // 3 = Between vertices 2 and 3
   // 4 = Between vertices 4 and 5
   // 5 = Between vertices 5 and 7
   // 6 = Between vertices 4 and 6
   // 7 = Between vertices 6 and 7
   // 8 = Between vertices 0 and 4
   // 9 = Between vertices 1 and 5
   // 10 = Between vertices 2 and 6
   // 11 = Between vertices 3 and 7
   // 
   // Pass the grid size for efficiency purposes.
   
   // Clear return values
   p1 = NULL;
   p2 = NULL;
   center = NULL;
   
   // Get the octant
   Octant* octant = getOctant( depth, grid_x, grid_y, grid_z );
   if ( octant == NULL ) {
      return;
   }
   
   // Set end points of the edge
   switch (edge)
   {
      case 0:
         p1 = octant->getVertex(0);
         p2 = octant->getVertex(1);                    
         break;
      case 1:
         p1 = octant->getVertex(1);
         p2 = octant->getVertex(3);
         break;
      case 2:
         p1 = octant->getVertex(0);
         p2 = octant->getVertex(2);
         break;
      case 3:
         p1 = octant->getVertex(2);
         p2 = octant->getVertex(3);
         break;
      case 4:
         p1 = octant->getVertex(4);
         p2 = octant->getVertex(5);
         break;
      case 5:
         p1 = octant->getVertex(5);
         p2 = octant->getVertex(7);
         break;
      case 6:
         p1 = octant->getVertex(4);
         p2 = octant->getVertex(6);
         break;
      case 7:
         p1 = octant->getVertex(6);
         p2 = octant->getVertex(7);
         break;
      case 8:
         p1 = octant->getVertex(0);
         p2 = octant->getVertex(4);
         break;
      case 9:
         p1 = octant->getVertex(1);
         p2 = octant->getVertex(5);
         break;
      case 10:
         p1 = octant->getVertex(2);
         p2 = octant->getVertex(6);
         break;
      case 11:
         p1 = octant->getVertex(3);
         p2 = octant->getVertex(7);
         break;
      default:
         break;
   };
   
   // Compute center point grid coordinates
   VectorObj center_grid_coords = (p1->getGridCoords() + p2->getGridCoords()) / 2.0;
   
   // Set center point, if it exists
   center = getVertex( center_grid_coords[0], center_grid_coords[1], center_grid_coords[2] );
}   

//-----------------------------------------------------
void
Octree::getOctantFace( 
   int depth, int grid_size, int grid_x, int grid_y, int grid_z, 
   int face, Vertex*& p1, Vertex*& p2, Vertex*& p3, Vertex*& p4, Vertex*& center )
{
  (void) grid_size;
   // Get a face of the given octant.
   // Returns the four vertices of the edge (where the odd/even vertices are diagonal
   // from each other) and the center point if it exists (NULL otherwise).
   //
   // Face index:
   // 0 = left, 1 = right, 2 = back, 3 = front, 4 = bottom, 5 = top

   // Get the octant
   Octant* octant = getOctant( depth, grid_x, grid_y, grid_z );

   // Get the vertex indices of the face
   VertexIndexVector vertex_indices;
   vertex_indices.resize(4);
   getVertexIndicesForFace( face, vertex_indices );
   
   // Set the vertices
   p1 = octant->getVertex( vertex_indices[0] );
   p2 = octant->getVertex( vertex_indices[1] );
   p3 = octant->getVertex( vertex_indices[2] );
   p4 = octant->getVertex( vertex_indices[3] );

   // Compute center grid coordinates
   VectorObj center_grid_coords = (p1->getGridCoords() + p3->getGridCoords()) / 2.0;

   // Set the center point, if it exists
   center = getVertex( center_grid_coords[0], center_grid_coords[1], center_grid_coords[2] );
}

//-----------------------------------------------------
void
Octree::getVertexIndicesForFace( int face, VertexIndexVector& vertex_indices )
{
   // Generically return the 4 vertex indices for the given face index:
   // Indices are returned in adjoining order (i.e. odd/even indices are
   // diagonal from each other).
   // 
   // Face index:
   // 0 = left, 1 = right, 2 = back, 3 = front, 4 = bottom, 5 = top
   //
   // Vertex index:
   // 0 = Corner (0, 0, 0)
   // 1 = Corner (1, 0, 0)
   // 2 = Corner (0, 0, 1)
   // 3 = Corner (1, 0, 1)
   // 4 = Corner (0, 1, 0)
   // 5 = Corner (1, 1, 0)
   // 6 = Corner (0, 1, 1)
   // 7 = Corner (1, 1, 1)
   // 8 = Center (0.5, 0.5, 0.5)
   
   switch (face)
   {
      case 0:
         vertex_indices[0] = 0;
         vertex_indices[1] = 2;
         vertex_indices[2] = 6;
         vertex_indices[3] = 4;
         break;
      case 1:
         vertex_indices[0] = 1;
         vertex_indices[1] = 3;
         vertex_indices[2] = 7;
         vertex_indices[3] = 5;
         break;
      case 2:
         vertex_indices[0] = 0;
         vertex_indices[1] = 1;
         vertex_indices[2] = 5;
         vertex_indices[3] = 4;
         break;
      case 3:
         vertex_indices[0] = 2;
         vertex_indices[1] = 3;
         vertex_indices[2] = 7;
         vertex_indices[3] = 6;
         break;
      case 4:
         vertex_indices[0] = 0;
         vertex_indices[1] = 1;
         vertex_indices[2] = 3;
         vertex_indices[3] = 2;
         break;
      case 5:
         vertex_indices[0] = 4;
         vertex_indices[1] = 5;
         vertex_indices[2] = 7;
         vertex_indices[3] = 6;
         break;
   };
}

//-----------------------------------------------------
void
Octree::getEdgeIndicesForFace( int face, EdgeIndexList& edge_indices )
{
   // Generically return the 4 edge indices for the given face index:
   // 
   // Face index:
   // 0 = left, 1 = right, 2 = back, 3 = front, 4 = bottom, 5 = top
   //
   // Edge index:
   // 0 = Between vertices 0 and 1
   // 1 = Between vertices 1 and 3
   // 2 = Between vertices 0 and 2
   // 3 = Between vertices 2 and 3
   // 4 = Between vertices 4 and 5
   // 5 = Between vertices 5 and 7
   // 6 = Between vertices 4 and 6
   // 7 = Between vertices 6 and 7
   // 8 = Between vertices 0 and 4
   // 9 = Between vertices 1 and 5
   // 10 = Between vertices 2 and 6
   // 11 = Between vertices 3 and 7
   
   // Clear return list
   edge_indices.clear();
   
   switch (face)
   {
      case 0:
         edge_indices.push_back( 2 );
         edge_indices.push_back( 6 );
         edge_indices.push_back( 8 );
         edge_indices.push_back( 10 );
         break;
      case 1:
         edge_indices.push_back( 1 );
         edge_indices.push_back( 5 );
         edge_indices.push_back( 9 );
         edge_indices.push_back( 11 );
         break;
      case 2:
         edge_indices.push_back( 0 );
         edge_indices.push_back( 4 );
         edge_indices.push_back( 8 );
         edge_indices.push_back( 9 );
         break;
      case 3:
         edge_indices.push_back( 3 );
         edge_indices.push_back( 7 );
         edge_indices.push_back( 10 );
         edge_indices.push_back( 11 );
         break;
      case 4:
         edge_indices.push_back( 0 );
         edge_indices.push_back( 1 );
         edge_indices.push_back( 2 );
         edge_indices.push_back( 3 );
         break;
      case 5:
         edge_indices.push_back( 4 );
         edge_indices.push_back( 5 );
         edge_indices.push_back( 6 );
         edge_indices.push_back( 7 );
         break;
      default:
         break;
   };
}

//-----------------------------------------------------
void
Octree::createEdge( Vertex* v1, Vertex* v2, edge_color_type edge_color )
{
   // Create an edge between the given two vertices, with the given color.
   
   // Check if the edge already exists
   if ( v1->getEdge( v2 ) != NULL ) {
      return;
   }
   
   // Create the edge
   Edge* new_edge = new Edge();
   
   // Set the vertices and color in the edge
   new_edge->setVertex1( v1 );
   new_edge->setVertex2( v2 );
   new_edge->setColor( edge_color );
   
   // Add this edge to the vertices
   v1->getEdgeList().push_back( new_edge );
   v2->getEdgeList().push_back( new_edge );
   
   // Add the edge to the octree
   m_edges.push_back( new_edge );
} 

//-----------------------------------------------------
void
Octree::clear()
{
   // Clear contents of this octree
      
   // Destroy the octants
   if ( m_depth > 0 )
   {
      for ( int i = m_depth; i >= 0; --i ) {
         OctantHashMap& hash_map = m_hash_map_vector[i];
         OctantHashMapIterator octant_iter = hash_map.begin();
         for ( ; octant_iter != hash_map.end(); ++octant_iter )
         {
            delete octant_iter->second;
         }
      }
      m_hash_map_vector.clear();
   }
   
   // Destroy the vertices
   VertexHashMapIterator vertex_iter = m_vertices.begin();
   for ( ; vertex_iter != m_vertices.end(); ++vertex_iter )
   {
      delete vertex_iter->second;
   }
   m_vertices.clear();
   
   // Destroy the edges
   EdgeVectorIterator edge_iter = m_edges.begin();
   for ( ; edge_iter != m_edges.end(); ++edge_iter )
   {
      delete *edge_iter;
   }
   m_edges.clear(); 

   m_depth = 0;
   m_leaf_width = 0.0;
   m_origin.clear();
}
