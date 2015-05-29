#ifndef TET_MESH_COLLISION_INTERFACE_H
#define TET_MESH_COLLISION_INTERFACE_H

#include <vector>
// comment out this line if you do not want use any acceleration structure for debuging purpose
//#define BRUTE_FORCE

class DeformBVHTree;
class TetCollider
{
public:
  /// \brief Float
  /// change double to float if the positions are of type float
  typedef double Float;
  /// \brief vert_tet_collision_threshold_
  /// Controls the threshold for deciding vertex-tet collision
  /// set to a positive small value like 0.01
  /// Larger value makes threshold distance for collision larger
  Float vert_tet_collision_threshold_;
  /// \brief vert_tet_collision_force_
  /// Controls how the amount of position correction to handle collision
  /// set to a positive small value like 0.01
  /// usually set to a value greater than vert_tet_collision_threshold_
  /// Larger value make correction position offset larger
  Float vert_tet_collision_force_;

  /// \brief edge_tri_collision_threshold_
  /// similar to vert-tet collision parameter
  Float edge_tri_collision_threshold_;
  Float edge_tri_collision_force_;
  /**
   * @brief TetCollider use this constructor if you are only doing vertex-tet collision
   * @param vert list of vertex positions
   * @param v_num vertex numver
   * @param tet list of vertex indices at each tetrahedron
   * @param tet_num number of tetrahedron
   */
  TetCollider(Float* vert, int v_num, int* tet, int tet_num);
  /**
   * @brief TetCollider use this constructor if you are doing both vert-tet and edge-tri colllision
   * @param vert list of vertex positions
   * @param v_num vertex numver
   * @param tet list of vertex indices at each tetrahedron
   * @param tet_num number of tetrahedron
   * @param edge  list of two verticies on edch edge
   * @param edge_num number of edges
   */
  TetCollider(Float* vert, int v_num, int* tet, int tet_num, int *edge, int edge_num);
  ~TetCollider();

  /**
   * @brief HandleCollision call this function at each time step to resolve collision
   * @param pos the current position of all verticies
   * @param new_pos the (possibly) collision-free positions after collision handling
   * @param kMaxIteration maximum number of iteration to handle collision
   */
  void HandleCollision(Float* pos, Float* new_pos, const int kMaxIteration = 100);

  //------------------------------------------------------------------------------
  // class internals you don't need to worry about if the code runs without problem
  void Init(Float* vert, int v_num, int* tet, int tet_num, int *edge = NULL, int edge_num = 0);
  Float (*vert_)[3];
  int v_num_;
  int (*tet_)[4];
  int tet_num_;
  int (*edge_)[2];
  int edge_num_;

  std::vector<Float> tmp_buffer_;
  std::vector<int> count_;
  DeformBVHTree* bvh_;
};

#endif // TET_MESH_COLLISION_INTERFACE_H
