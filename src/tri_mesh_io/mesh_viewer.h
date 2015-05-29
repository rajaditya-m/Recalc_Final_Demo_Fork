#ifndef MESHVIEWER_H
#define MESHVIEWER_H
#include <vector>
#include "input_handler.h"
#include "global.h"
#include "timer.h"

struct MeshFileType
{
  enum {
    kTetGen,
    kVega,
    kObj,
    kOff,
  };
};

class TetrahedralMesh;
class TriangularMesh;
class TetGenMeshIO;
class VegaTetMeshIO;
class ObjMeshIO;
class OffMeshIO;

template <int type>
struct MeshChooser {
};

template <>
struct MeshChooser<MeshFileType::kTetGen> {
  typedef TetrahedralMesh Mesh;
  typedef TetGenMeshIO MeshIO;
};

template <>
struct MeshChooser<MeshFileType::kVega> {
  typedef TetrahedralMesh Mesh;
  typedef VegaTetMeshIO MeshIO;
};

template <>
struct MeshChooser<MeshFileType::kObj> {
  typedef TriangularMesh Mesh;
  typedef ObjMeshIO MeshIO;
};

template <>
struct MeshChooser<MeshFileType::kOff> {
  typedef TriangularMesh Mesh;
  typedef OffMeshIO MeshIO;
};

class QTimer;
template <int file_type>
class MeshViewer : public InputHandler
{
public:
  typedef typename MeshChooser<file_type>::Mesh Mesh;
  typedef typename MeshChooser<file_type>::MeshIO MeshIO;
  virtual int HandleKeyPress(QKeyEvent* e);
  virtual int Render();
  virtual int Idle();
  MeshViewer(Mesh* mesh, std::vector<std::string>& mesh_files);
  virtual ~MeshViewer() {}
  void DumpRenderResultToPng(const char* folder = NULL);
  void set_frame_rate(double frame_rate);
  void set_output_path(std::string path);
  void LoadMesh(std::string mesh_file);
  bool play_animation_;
  Mesh* mesh_;
  std::string output_path_;
  uint64_t prev_time_;
  bool dumping_render_result_;
  double frame_rate_;
  int current_mesh_idx_;
  std::vector<std::string> mesh_files_;
};

extern template class MeshViewer<MeshFileType::kTetGen>;
extern template class MeshViewer<MeshFileType::kVega>;
extern template class MeshViewer<MeshFileType::kObj>;
extern template class MeshViewer<MeshFileType::kOff>;

#endif // MESHVIEWER_H
