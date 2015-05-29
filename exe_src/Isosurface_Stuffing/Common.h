
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

#ifndef COMMON_H
#define COMMON_H

#include <list>
#include <vector>
//#include <hash_map>
#include <unordered_map>

/**
 **  This header file contains all the commonly used typedefs, defines, etc.
 **/

class Vertex;
class Edge;
class Octant;
class Tetrahedron;
class VectorObj;

// Enumerations
enum edge_color_type
{
   BLACK_EDGE = 0,
   RED_EDGE,
   BLUE_EDGE
};

enum tetrahedron_type
{
   BCC_TETRAHEDRON = 0,
   BISECTED_BCC_TETRAHEDRON,
   QUADRISECTED_BCC_TETRAHEDRON,
   HALF_PYRAMID_TETRAHEDRON
};

// Typedefs
typedef std::vector<Vertex*> VertexVector;
typedef VertexVector::iterator VertexVectorIterator;
typedef VertexVector::const_iterator VertexVectorConstIterator;

typedef std::list<Vertex*> VertexList;
typedef VertexList::iterator VertexListIterator;
typedef VertexList::const_iterator VertexListConstIterator;

typedef std::vector<Edge*> EdgeVector;
typedef EdgeVector::iterator EdgeVectorIterator;
typedef EdgeVector::const_iterator EdgeVectorConstIterator;

typedef std::list<Edge*> EdgeList;
typedef EdgeList::iterator EdgeListIterator;
typedef EdgeList::const_iterator EdgeListConstIterator;

typedef std::vector<Tetrahedron*> TetrahedronVector;
typedef TetrahedronVector::iterator TetrahedronVectorIterator;
typedef TetrahedronVector::const_iterator TetrahedronVectorConstIterator;

typedef std::vector<Octant*> OctantVector;
typedef OctantVector::iterator OctantVectorIterator;
typedef OctantVector::const_iterator OctantVectorConstIterator;

typedef std::pair<long, Octant*> OctantPair;
//typedef stdext::hash_map<long, Octant*> OctantHashMap;
typedef std::unordered_map<long, Octant*> OctantHashMap;
typedef OctantHashMap::iterator OctantHashMapIterator;
typedef OctantHashMap::const_iterator OctantHashMapConstIterator;

typedef std::vector<OctantHashMap> OctantHashMapVector;
typedef OctantHashMapVector::iterator OctantHashMapVectorIterator;
typedef OctantHashMapVector::const_iterator OctantHashMapVectorConstIterator;

typedef std::pair<long, Vertex*> VertexPair;
typedef std::unordered_map<long, Vertex*> VertexHashMap;
typedef VertexHashMap::iterator VertexHashMapIterator;
typedef VertexHashMap::const_iterator VertexHashMapConstIterator;

typedef std::list<int> EdgeIndexList;
typedef EdgeIndexList::iterator EdgeIndexListIterator;
typedef EdgeIndexList::const_iterator EdgeIndexListConstIterator;

typedef std::vector<int> VertexIndexVector;
typedef VertexIndexVector::iterator VertexIndexVectorIterator;
typedef VertexIndexVector::const_iterator VertexIndexVectorConstIterator;

typedef std::vector<VectorObj> GridCoordVector;
typedef GridCoordVector::iterator GridCoordVectorIterator;
typedef GridCoordVector::const_iterator GridCoordVectorConstIterator;

#endif
