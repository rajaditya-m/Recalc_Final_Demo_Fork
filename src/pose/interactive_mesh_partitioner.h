#ifndef INTERACTIVE_MESH_PARTITIONER_H
#define INTERACTIVE_MESH_PARTITIONER_H
#include <vector>
#include "input_handler.h"
#include "vector_lib.h"
//#include <vector>
#include "tet.h"

class Tet;
class InteractiveMeshPartitioner : public InputHandler
{
public:
  typedef dj::Vec3d Vec3;
  InteractiveMeshPartitioner(Tet* tet);
  virtual ~InteractiveMeshPartitioner() {}
  virtual int HandleMouseRelease(QMouseEvent *e);
  virtual int HandleKeyPress(QKeyEvent *e);
//  virtual int HandleKeyRelease(QKeyEvent *e) { return -1;}
//  virtual int HandleMouseMove(QMouseEvent*e) {return -1;}
  virtual int HandleMousePress(QMouseEvent *e);
  virtual int Render();
  enum {
    kNum = 3,
  };

  enum RenderMode {
    kRenderSelectedVertex = 1 << 0,
    kRenderPartitionedVertex = 1 << 1,
    kRenderUnSelectedVertex = 1 << 2,
    kRenderTetMesh = 1 << 3,
  };

  static const int kRenderModes[];
  static const int kRenderModeNum;

  enum {
    kDefault,
    kFirstLeftButton,
    kSecondLeftButton,
    kThirdLeftButton,
    kStateNum
  };
  enum VertexState {
    kUnSelected = -2,
    kSelected = -1,
    kPartitioned,
  };
  std::vector<int> vertex_group_;
  void SavePartition(const char* file_name);
  Vec3 ray_start_[kNum];
  Vec3 ray_end_[kNum];
  Vec3 clicked_points_[kNum];
  int current_group_num_;
  int render_mode_;
  int mouse_state_;
  Tet* tet_;
};

#endif // INTERACTIVE_MESH_PARTITIONER_H
