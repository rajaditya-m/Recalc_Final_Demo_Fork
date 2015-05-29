#ifndef	__WHMIN_TET_H__
#define __WHMIN_TET_H__
//#include <armadillo>
#include <vector>
#include <unordered_map>
#include "selectable_object.h"
#include "buffer_manager.h"
#include "MY_MATH.h"
//#include "VECTOR_TOOLS.h"

class QMouseEvent;
class TetSampler;
template<class Float> class ConjugateGradientSolver;
template<class Real> class TriangleMeshRender;
template <class Real> class TexturedTriangleMeshRenderer;
//bool HandleMouseMove(QMouseEvent * event);
//bool HandleMouseRelease(QMouseEvent * event);
//bool HandleMousePress(QMouseEvent * event, bool moving_body);
//void HandleKeyPress(unsigned char key);

namespace dj {
template <class T, int d, bool> class Vec;
//typedef Vec<double, 3, false> Vec3f;
typedef Vec<double, 3, false> Vec3d;
}

template <class T> class AffineTransformer;
class Skeleton;
class TetMeshSimulatorBridge;
class Tet : public SelectableObject<double> {
protected:
  BufferManager buffer_;

public:
  enum MeshRenderMode {
    kDefaultRendering = 0,
    kWireFrame = 1 << 1,
    kSmoothLighting = 1 << 2,
    kFlatLighting = 1 << 3,
    kSurfaceEdge = 1 << 4,
    kTransparent = 1 << 5,
    kParitionedMeshVertex = 1 << 6,
    kTets = 1 << 7,
    kNoRendering = 1 << 8,
  };
  static const int kRenderModeNum;
  static const int kRenderMode[];

  int render_mode_;

  int one_indexed_;
  Skeleton* skeleton_;
  std::vector<int> surface_vertices_;
  ConjugateGradientSolver<double>* cg_solver_;
  int selected_vertex_;
  //Vertex
  double ui_force_scaling_factor_;
  double internal_force_scaling_factor_;
  int max_incident_tet_;
  int* has_detail_;
  int vertex_num_;
  double* stiffness_;
  double* rest_pos_;
  double* X;
  double* tmp_X;
  double* render_X;
  double* strain_limiting_offset_;
  int* x_tag;
  int* closest_bone_;
  double* initial_offset_;
  double* old_X;
  double* V;
  double* temp;
  double* mass_;
  int* weight;
  double* velocity_;
   int domain_offset_toggle_;
  TriangleMeshRender<double>* surface_renderer_;
  std::vector<TexturedTriangleMeshRenderer<double>> textured_surface_renderer_;
  std::vector<double> surface_vert_normal_;
  std::vector<double> surface_vert_pos_;
  std::vector<double> surface_vert_texture_coord_;
  std::vector<int> surface_triangle_indices_; // surface triangle index in surface vertex only array
  std::vector<double> surface_vert_colors_;

  std::vector<std::vector<int> > adjacent_vertex_;
  std::vector<std::vector<int> > incident_edge_;
  std::vector<std::vector<int> > incident_triangle_;
  std::vector<std::vector<int> > incident_triangle_renderer_;
  std::vector<std::vector<int> > incident_triangle_renderer__2;
  std::vector<std::vector<int> > special_triangles_;
  std::vector<std::vector<int> > incident_tet_;
  std::vector<std::vector<int> > rank_in_tet_;
  double* tmp_vertex_pos_;

  int*	tet_;
  int tet_number;

  double* Dm;
  double* inv_Dm;
  double* Area_Coeffs;
  double* volume_;

  // triangles for rendering purpose
  int triangle_num_;
  double* triangle_vector_;
  double* triangle_bounding_box_;
  int* T;
  int* T_Stitch_;
  int last_tri;
  int* attached_tet_of_surface_triangle_;
  double* VN; //Vertex Normal
  double* TN; //Triangle Normal
  double* triangle_area_;
  int cubature_flip_counter_;

  int* edges_; //edge list
  double* rest_length_;
  int edge_num_;
  double* tet_weights_; // 1 / max(incident tet of verticies);
  double density_;
  bool position_changed_;
  std::vector<double> edge_stiffness_;

  TetMeshSimulatorBridge* inv_fem_;
  std::vector<int> surface_tets_;
  std::vector<int> is_boundary_tri_on_tet_;

  double tet_limit_threshold;
  TetSampler* higher_sampler_;
  TetSampler* lower_sampler_;
  double per_vertex_weight_;

  //  std::vector<int> is_on_skeleton_;
  void GetJointPositions(double* joint_pos);
  std::vector<int> attached_joint_of_vertex_;
  std::vector<int> constrainted_vert_vel_;
  std::vector<int> constrainted_vertex_;
  std::vector<int> is_constrainted_;
  std::vector<std::pair<int, double> > joints_;
  std::vector<double> joint_pos_;
  std::vector<double> joint_vel_;
  std::vector<int> verts_binded_to_skeleton_;
  std::vector<int> attached_tet_of_joints_;
  std::vector<double> barycentric_of_joint_;

  // Some more rendering data structure
  std::vector<int> interfaceTriangles_;
  std::vector<int> interface_vertices_global_;
  int surfaceVertOldSize_;

  std::vector<std::pair<int,int> > interface_vertices_;
  int interface_vertex_split_;

  Tet(const char *filename, double _limit_threshold, AffineTransformer<double>* affine_transform = NULL,  bool initialize_fem_module = true,bool load_interface = false, int split=0);
  virtual ~Tet();

  void setCubatureFlipCounter(int i) { cubature_flip_counter_ = i;  }

  virtual int Select(double* ray_start, double* ray_end, double* clicked_world_pos, double*selected_pos);
  virtual bool GetVertexPosition(double* pos);

  int GetUIForce(double* force);

  virtual void UpdateRenderData();

  void UpdateVegaVertexOffset();
  void NextRenderMode() {
    render_mode_ = (render_mode_ + 1) % kRenderModeNum;
  }
  void PrevRenderMode() {
    render_mode_--;
    if (render_mode_ < 0) render_mode_ = kRenderModeNum - 1;
  }
  void MoveSkeletonJoint(int joint_idx, int angle_idx, double angle_change_in_radian);

  void SaveConstrainedVertex(const char* file_name);

  void SaveBoneVertex(const char* file_name);
  void SaveMeshPartition(const char* file_name);
  void LoadMeshPartition(const char* file_name);
  virtual void EnalbeSurfaceRenderingWithVBO(const char* shader_source);
  virtual void EnalbeTexturedSurfaceRendering(const char* shader_source);


  void Reset();

  int VertexNum() {
    return vertex_num_;
  }
  int TetNum() {
    return tet_number;
  }
  double* VertexArray() {
    return X;
  }
  int* TetIndex() {
    return tet_;
  }


  // assume normal is up-to-date
  virtual void UpdateSurfaceMeshVBO();
  virtual void UpdateTexturedSurfaceMeshVBO();
  void MassSpringSimulation(double dt);
  void AllocateVertexData(int v_num);
  void AllocateTriangleData(int tri_num);
  void AllocateTetData(int tet_num);
  void AllocateEdgeData(int e_num);
  void Simulate(double dt);

  void ComputeSurfaceTriangleArea();

  void ExportTetMesh(const char* file_name);

  void ImportTetMesh(const char* file_name);

  void set_frame_num(int new_frame_num);
  inline void set_current_frame(int new_current_frame);
  inline int current_frame();
  void InitializeMass();

  //------------------------------------------------------------------------------
  // Topology and geomtry Initialization functions
  //------------------------------------------------------------------------------
  std::vector<int> GetIncidentTetStat();

  void FindBodyAttachedTetOfJoints(Skeleton *skeleton);
  void AttachSkeleton(Skeleton *skeleton, bool find_attached_tet = true);
  void AttachSkeleton(Skeleton *skeleton, Skeleton* rest_skeleton, const char *file_name);
  void ReinitializeOffsetToSkeleton();
  void Build_Edge_From_Tets();
  void Quick_Sort_RE( int a[], int l, int r);
  int Quick_Sort_Partition_RE( int a[], int l, int r);
  void ComputeRestDeformationGradientAndVolume();
  void Build_Boundary_Triangles();
  void QuickSort( int a[], int l, int r);
  int QuickSort_Partition( int a[], int l, int r);
  void Build_TN();
  void Build_VN();
  double ComputeAvgEdgeLength();

  //------------------------------------------------------------------------------
  // IO functions
  //------------------------------------------------------------------------------
  void Write_File(const char *file_name);
  void Read_File(const char *file_name);
  void Read_Original_File(const char *prefix_name);
  void Write_Original_File(const char *name);
  void Load_SMF_File(char *name);
  void Output_Mesh(char *file_name);
  void Input_Mesh(char *file_name);
  void Load_PLY_File(char *name);
  void Output_POLY_File(const char *name);
  void Save_OBJ_File(const char *name);
  void Save_PLY_File(char *name);
  void Output_STL_File(char *name);
  void SavePosition(const char* file_name);
  void LoadPosition(const char* file_name);
  void LoadPosition(std::vector<double> &pos);

  //------------------------------------------------------------------------------
  // Rendering functions
  //------------------------------------------------------------------------------
  void Render(int render_mode = kDefaultRendering, double *pos = NULL);

  //------------------------------------------------------------------------------
  // Strain limiting functions
  //------------------------------------------------------------------------------
  void ReinitializeSkinningWeight() {
    //    for (int v = 0; v < vertex_num_; ++v) {
    //      SkeletonNode* closest_node = skeleton_->GetNode(closes_joint_[v]);
    //      dj::Mat3f parent_rotation(closest_node->parent()->rotation_matrix());
    //      dj::Mat3f rotation_matrix(closest_node->rotation_angle_matrix_);
    //      rotation_matrix = parent_rotation * rotation_matrix;
    //      rotation_matrix.Transpose();
    //      Vec3f offset = Vec3f(X + v * 3) - Vec3f(closest_node->parent()->world_pos());
    //      dj::FixedVecWrapper<double, 3> initial_offset(initial_offset_ + v * 3);
    //      initial_offset = rotation_matrix * offset;
    //    }
  }
  // User interaction and skeleton related functions
  //------------------------------------------------------------------------------
  void ApplySkeletonTransformationToVertex();
  void GetConstrainVertexVelocity();
  int SelectVertex(const double* start, const double* end, double *select_vertex_pos, double* offset);

  void ComputeTriVecNormalAndBoundingBox(double threshold);
private:
  void BuildIncidentTet(void);
};








#endif  //__WHMIN_TET_H__
