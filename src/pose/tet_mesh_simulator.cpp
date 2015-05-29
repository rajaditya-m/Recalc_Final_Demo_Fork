#include <stdio.h>
#include <fstream>
//#include "shape_op_tet.h"
#include "tet_mesh_simulator.h"
#include "global.h"
#include "qt_object_selector.h"
#include "skeleton.h"
#include "skeleton_node.h"
#include "affine_transformer.h"
#include "tet.h"
#include "open_gl_qt.h"
#include "config_file.h"
#include "tet_gen_mesh_io.h"

namespace pose_simulator {
int render_mode = 0;
//ShapeOpTet* shape_op_tet = NULL;

using namespace global;
void CreateBody() {
  L("Start creating bodies");
  std::vector<std::string>* transform_cmd = conf.Get<std::vector<std::string>* >("affine_transform");
//  transform_cmd->push_back("range 1");
//  transform_cmd->push_back("centerize");
//  transform_cmd->push_back("min_y 0");
  AffineTransformer<double> transform(*transform_cmd);
//  const char* mesh_file = DATA_DIRECTORY "tetgen_mesh/P";
//  const char* mesh_file = DATA_DIRECTORY "tetgen_mesh/cube10k";
//  const char* mesh_file = DATA_DIRECTORY "tetgen_mesh/hand";
//  const char* mesh_file = DATA_DIRECTORY "beam_8x4x20/beam_8x4x20";
//  shape_op_tet = new ShapeOpTet(mesh_file, &transform, false);
//  global::current_body = shape_op_tet;
  global::current_body = new Tet(GetMeshFile(), 0, &transform);
//  global::current_body->Save_OBJ_File("/Users/dj/octopus.obj");
//  exit(0);
#if 1
  // TODO delete
#else
  bodies.back()->AttachSkeleton(skeleton);
#endif
  L("finish creating body");
}

void CreateSkeleton() {
  L("Create Skeleton");
  //  P(global::skeleton_file);
  skeleton = new Skeleton(global::skeleton_file.c_str(), true);
}


void Init() {
//  CreateSkeleton();
  CreateBody();
  global::gl->InstallHandler(new QtObjectSelector<double>(global::current_body));
  //  exit(0);
  //  // TODO disable
  if (0) {
    std::ifstream in(DATA_DIRECTORY "octopus/vert_partition.txt");
    ASSERT(in.is_open());
    int v_num;
    in >> v_num;
    ASSERT(v_num == global::current_body->vertex_num_);
    for (int v = 0; v < v_num; ++v) {
      in >> global::current_body->is_constrainted_[v];
    }
    in.close();
  }
}

} // namespace pose_simulator

using namespace pose_simulator;

TetMeshSimulator::TetMeshSimulator()
  : InputHandler(2) {
  Init();
}

TetMeshSimulator::~TetMeshSimulator() {
}

int TetMeshSimulator::Idle() {
  global::current_body->Simulate(global::time_step);
//  shape_op_tet->Simulate(global::time_step);
  bool export_image = conf.Get<int>("export image");
  bool export_obj = conf.Get<int>("export obj");
  if (export_image || export_obj) {
    static int count = 0;
    global::gl->updateGL();
    char file[1 << 10];
    int step_per_frame = int(1.0 / global::time_step / 30.0);
    if (step_per_frame == 0) step_per_frame = 1;
    if (count % step_per_frame == 0) {
      int frame_id = count / step_per_frame;
#if defined(__APPLE__)
      const char folder[] = "/Volumes/ram";
#elif defined(_WIN32) || defined(_WIN64)
      const char folder[] = "C:/tmp/render";
#else
      const char folder[] = "/tmp/log";
#endif
      if (export_image) {
        sprintf(file, "%s/%05d.png", folder, frame_id);
        global::gl->ScreenShot(file, "png");
      }
      if (export_obj) {
        sprintf(file, "%s/obj/%05d.obj", folder, frame_id);
        global::current_body->Save_OBJ_File(file);
      }
    }
    count++;
  }
  return kNotHandled;
}

int TetMeshSimulator::Render() {
  global::current_body->Render(Tet::kRenderMode[render_mode], current_body->X);
  //  global::skeleton->Render();
  return kNotHandled;
}

int TetMeshSimulator::HandleKeyPress(QKeyEvent *e) {
  switch (e->key()) {
    case Qt::Key_R:
      render_mode = (render_mode + 1) % Tet::kRenderModeNum;
      break;
    case Qt::Key_X: {
      if (0) {
        std::vector<double> verts(global::current_body->X, global::current_body->X + global::current_body->vertex_num_ * 3);
        std::vector<int> tets(global::current_body->tet_, global::current_body->tet_ + global::current_body->tet_number * 4);
        TetGenMeshIO::Instance()->Write("/Users/dj/octopus", verts, tets);
        L("saved");
        break;
      }
      char file_name[512];
      sprintf(file_name, "%s/body%d.obj", DATA_DIRECTORY, 0);
      L("Save body " + std::string(file_name));
      global::current_body->Save_OBJ_File(file_name);
      break;
    }
    default:
      return -1;
  }
  return kPriority_;
}

