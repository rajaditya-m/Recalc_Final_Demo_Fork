#include <fstream>
#include <Eigen/Dense>
#include "mixed_multi_domain_simulator.h"
#include "rainbow_color.h"
#include "mixed_multi_domain_tet.h"
#include "multi_domain_tet.h"
#include "skeleton_node.h"
#include "skeleton.h"
#include "camera.h"
#include "global.h"
#include "open_gl_qt.h"
#include "affine_transformer.h"
#include "tet_mesh_simulator.h"
#include "affine_transformer.h"
#include "string_formatter.h"
#include "objMesh.h"
#include "objMeshRender.h"
#include "qt_object_selector.h"
#include "config_file.h"

//#undef DATA_DIRECTORY
//#define DATA_DIRECTORY "/tmp/log/"
namespace mixed_multi_domain_simulator {
enum UIMode {
  kDefault,
  kCut
};
enum RenderMode {
  kEdge = 0,
  kSurface = 1,
};
//Eigen::Vector3d pa(-2.74963, 2.19545, 1.82949);
//Eigen::Vector3d pb(0.192739, 0.298575, 0.991291);
//Eigen::Vector3d pc(-0.0690679, 0.0892717, 0.750821);

//Eigen::Vector3d pa(-0.781474, 0.429062, 0.340563);
//Eigen::Vector3d pb(9785.43, -1274.3, -2061.64);
//Eigen::Vector3d pc(8825.78, -3933.79, -3473.16);

Eigen::Vector3d pa(-3.33686, 1.45635, 1.20643);
Eigen::Vector3d pb(9857.42, -2876.28, -510.322);
Eigen::Vector3d pc(9404.48, -4400.79, -499.782);


int render_mode = kSurface;
bool render_cut_plane = false;
MixedMultiDomainTet* multi_tet;

using namespace global;

void CreateBody() {
  L("Start creating bodies");
  std::vector<std::string> cmd;
  //  cmd.push_back("translate 0 0.5 0");
  cmd.push_back("centerize");
  cmd.push_back("range 1");
  cmd.push_back("min_x 0.0");
  cmd.push_back("min_y 0.0");
  cmd.push_back("min_z 0.0");
  AffineTransformer<double> transform(cmd);
  //      multi_tet = new MixedMultiDomainTet(DATA_DIRECTORY "tetgen_mesh/cube", &transform);
  //  multi_tet = new MixedMultiDomainTet(DATA_DIRECTORY "tetgen_mesh/beam_2x3x2", &transform);
  //  multi_tet = new MixedMultiDomainTet(DATA_DIRECTORY "tetgen_mesh/beam_2x2x2", &transform);
  //    multi_tet = new MixedMultiDomainTet(DATA_DIRECTORY "tetgen_mesh/P", &transform);
  multi_tet = new MixedMultiDomainTet(GetMeshFile());
  //          multi_tet = new MixedMultiDomainTet(DATA_DIRECTORY "tetgen_mesh/one_tet");
  //      multi_tet = new MixedMultiDomainTet(DATA_DIRECTORY "tetgen_mesh/two_tet");
  //  std::vector<std::string>* transform_cmd = pose_conf.Get<std::vector<std::string>* >("affine_transform");
  // todo uncomment
  //  multi_tet->AttachSkeleton(global::skeleton);
  L("finish creating body");
}

void CreateSkeleton() {
  return;
  L("Create Skeleton");
  skeleton = new Skeleton(global::skeleton_file.c_str(), true);
}

void Test() {
  if (0) {
    multi_tet->Cut(pa, pb, pc);
  }
  if (0) {
    Eigen::Vector3d a(-1.08, 0.55, 0.15);
    Eigen::Vector3d b(1.08, 0.55, 0.15);
    Eigen::Vector3d dir(0.0, -0.45, 0);
    multi_tet->Cut(a, b, dir);
  }
  if (0) {
    Eigen::Vector3d a(-1.08, 0.55, 0.25);
    Eigen::Vector3d b(1.08, 0.55, 0.25);
    Eigen::Vector3d dir(0.0, -0.45, 0);
    multi_tet->Cut(a, b, dir);
  }

  if (0) {
    Eigen::Vector3d a(-1.00, -0.03, -1.15);
    Eigen::Vector3d b(-1.00, -0.03, +1.15);
    Eigen::Vector3d dir(3.0, 0, 0);
    multi_tet->Cut(a, b, dir);
  }

  if (0) {
    Eigen::Vector3d a(-1.08, 1.55, 0.15);
    Eigen::Vector3d b(1.08, 1.55, 0.15);
    Eigen::Vector3d dir(0.0, -5, 0);
    multi_tet->Cut(a, b, dir);
  }
  if (0) {
    Eigen::Vector3d a(-1.08, 1.25, 0.25);
    Eigen::Vector3d b(1.08, 1.25, 0.25);
    Eigen::Vector3d dir(0.0, -0.75, 0);
    multi_tet->Cut(a, b, dir);
  }
  if (0) {
    Eigen::Vector3d a(0.8, 1.2, -1);
    Eigen::Vector3d b(0.8, 1.2, 1);
    Eigen::Vector3d dir(0.0, -0.8, 0);
    multi_tet->Cut(a, b, dir);
  }

  if (0) {
    Eigen::Vector3d a(-0.5, 1, -4);
    Eigen::Vector3d b(1.5, -1, -4);
    Eigen::Vector3d dir(0.0, 0, 6);
    multi_tet->Cut(a, b, dir);
  }
}

void Init() {
  CreateSkeleton();
  CreateBody();
  global::current_body = multi_tet;
  global::gl->InstallHandler(new QtObjectSelector<double>(global::current_body));
  auto GetSubspaceFileName = [](int part_id) -> const char* {
    static std::string name;
    name = dj::Format("%s/modal_basis/partition_%d.basis.bin", GetDataFolder(), part_id);
    return name.c_str();
  };
  if (1) {
    multi_tet->LoadPartitionInfo(GetPartitionFolder());
    multi_tet->LoadSubspace(GetSubspaceFileName, MultiDomainTet::kBinary);
    //    multi_tet->AssembleGlobalBasis();
    multi_tet->LoadCubature(dj::Format("%s/cubature", GetDataFolder()).c_str());
    std::vector<int>* full_domains = conf.Get<std::vector<int>*>("full simulation domains");
    std::unordered_set<int> full_sim_domain(full_domains->begin(), full_domains->end());
    multi_tet->PrecomputeFastSandwichTransform();
    multi_tet->SetFullSimulationDomains(full_sim_domain);
  }
  std::vector<int> fixed_domain = *conf.Get<std::vector<int>*>("fixed domain");
  if (fixed_domain[0] != -1) {
    multi_tet->SetFixedDomains(std::set<int>(fixed_domain.begin(), fixed_domain.end()));
  }

  //  multi_tet->GenerateCubicaStylePartition(DATA_DIRECTORY "armadillo/cubica_partition"); L("cubica partition generated"); exit(0);
  //  multi_tet->SetConstrainedVertex();
  //    multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/ascMat");
  //    multi_tet->SaveConstrainedVertex(DATA_DIRECTORY "armadillo/armadillo_150k.fixed_vertex.txt");
  Test();
  multi_tet->Build_TN();

  if (conf.Get<int>("embed obj mesh")) {
    std::string obj_file = dj::Format("%z/../../3D Model of Octopus/OCTOPUS/OCTOPUS/octopus.3.obj", DATA_DIRECTORY);
    static ObjMesh obj_(obj_file.c_str(), ObjMesh::ASCII, false);
    static ObjMeshRender obj_render(&obj_);
    int textureMode = OBJMESHRENDER_GL_USEANISOTROPICFILTERING | OBJMESHRENDER_GL_USEMIPMAP | OBJMESHRENDER_LIGHTINGMODULATIONBIT;
    obj_render.loadTextures(textureMode);
    //    multi_tet->UpdateObjMeshPosition();
    std::string mapping_folder = dj::Format("%z/obj_mapping", GetModelFolder());
    multi_tet->EmbededObjMesh(&obj_, &obj_render, true, mapping_folder.c_str());
    multi_tet->EnalbeTexturedSurfaceRendering(DATA_DIRECTORY "../shader/octopus");
  } else {
    multi_tet->EnalbeSurfaceRenderingWithVBO(DATA_DIRECTORY "../shader/phong");
  }
  if (conf.Get<int>("export video")) {
    global::gl->EnableVideoRecording();
  }
  //  multi_tet->EnalbeTexturedSurfaceRendering(DATA_DIRECTORY "../shader/phong");
}

} // namespace pose_simulator

using namespace mixed_multi_domain_simulator;


int MixedMultiDomainSimulator::HandleMousePress(QMouseEvent *e) {
  if ((e->modifiers() & Qt::ControlModifier) && e->button() == Qt::RightButton) {
    QPoint pos = e->pos();
    //    double start[3], end[3];
    //    double clicked_world_pos[3];
    //    GetPixelWorldPosition(pos.x(), pos.y(), clicked_world_pos);
    GetPixelWorldPosition(pos.x(), pos.y(), clicked_pos_);
    depth_ = GetPixelDepth(pos.x(), pos.y());
    //    GetSelectionRay(pos.x(), pos.y(), start, end);
    //    int intersect = multi_tet->Intersect(start, end, clicked_world_pos, clicked_pos_);
    //    if (intersect >= 0) {
    ui_mode_ = kCut;
    //    }
  }
  return kNotHandled;
}


int MixedMultiDomainSimulator::HandleMouseRelease(QMouseEvent *e) {
  if (ui_mode_ == kCut) {
    QPoint pos = e->pos();
    MultiDomainTet::Vec3 second_click_pos;
    GetPixelWorldPosition(pos.x(), pos.y(), depth_, &second_click_pos[0]);
    //    GetSelectionRay(pos.x(), pos.y(), &start[0], &end[0]);
    //    MultiDomainTet::Vec3 dir = MultiDomainTet::MapVec3(clicked_pos_) - start;
    pa[0] = gl->camera()->eye_pos()[0];
    pa[1] = gl->camera()->eye_pos()[1];
    pa[2] = gl->camera()->eye_pos()[2];
    pb = MultiDomainTet::Vec3(clicked_pos_);
    pc = second_click_pos;
//    PVEC(pa);
//    PVEC(pb);
//    PVEC(pc);
    multi_tet->Cut(pa, pb, pc);
//    L("cut");
    ui_mode_ = kDefault;
    gl->updateGL();
  }
  return kNotHandled;
}

int MixedMultiDomainSimulator::HandleKeyPress(QKeyEvent *e) {
  switch (e->key()) {
    case Qt::Key_E:
      render_cut_plane = !render_cut_plane;
      break;
    case Qt::Key_R:
      render_mode = 1 - render_mode;
      break;
    default:
      break;
  }
  return kNotHandled;
}

int MixedMultiDomainSimulator::Render() {
  profiler.Start("render");
  if (render_mode == kSurface) {
    multi_tet->RenderSurface();
  } else {
    multi_tet->Render();
  }
  if (render_cut_plane) {
    glEnable(GL_LIGHTING);
    glColor3fv(kBlue());
    glBegin(GL_TRIANGLES);
    Vertex3v(&pa[0]);
    Vertex3v(&pb[0]);
    Vertex3v(&pc[0]);
    glEnd();
  }
  profiler.End("render");
  return kNotHandled;
}

int MixedMultiDomainSimulator::Idle() {
  static double total_simulation_time = conf.Get<double>("total simulation time");
  static unsigned int steps = 0;
  for (int i = 0; i < global::simulation_step_per_idle; ++i) {
    if (total_simulation_time > 0 && steps * global::time_step > total_simulation_time) {
      global::simulate = false;
      //      exit(0);
      continue;
    }
    steps++;
    multi_tet->MixSimulaitonOneStepWithCubature(global::time_step);
    bool export_image = conf.Get<int>("export image");
    bool export_obj = conf.Get<int>("export obj");
    bool export_video = conf.Get<int>("export video");
    if (export_image || export_obj || export_video) {
      static int count = 0;
      char file[1 << 10];
      int step_per_frame = int(1.0 / global::time_step / 30.0);
      if (step_per_frame == 0) step_per_frame = 1;
      if (count % step_per_frame == 0) {
        if (export_video) {
          global::gl->RecordNextFrame();
        }
        int frame_id = count / step_per_frame;
#if defined(__APPLE__)
        const char folder[] = "/Volumes/ram";
        //        const char folder[] = "/tmp/png";
#elif defined(_WIN32) || defined(_WIN64)
        const char folder[] = "C:/tmp/render";
#else
        const char folder[] = "/tmp/log";
#endif
        if (export_image) {
          if (!export_video) global::gl->updateGL();
          sprintf(file, "%s/%05d.png", folder, frame_id);
          global::gl->ScreenShot(file, "png");
        }
        if (export_obj) {
          sprintf(file, "%s/obj/%05d.obj", folder, frame_id);
          multi_tet->Save_OBJ_File(file);
        }
      }
      count++;
    }
  }
  //  multi_tet->UpdateTexturedSurfaceMeshVBO();
  return kNotHandled;
}

MixedMultiDomainSimulator::MixedMultiDomainSimulator() {
  Init();
  ui_mode_ = kDefault;
}

MixedMultiDomainSimulator::~MixedMultiDomainSimulator() {

}
