#include <fstream>
#include <unordered_set>
#include "multi_domain_subspace_simulator.h"
#include "multi_domain_tet.h"
#include "skeleton_node.h"
#include "skeleton.h"
#include "global.h"
#include "open_gl_qt.h"
#include "tet_mesh_simulator.h"
#include "affine_transformer.h"
#include "string_formatter.h"
#include "config_file.h"
#include "obj_mesh_io.h"

#include "objMesh.h"
#include "objMeshRender.h"
//#undef DATA_DIRECTORY
//#define DATA_DIRECTORY "/tmp/log/"
namespace multi_domain_simulator {
MultiDomainTet* multi_tet;
int render_mode = 0;

using namespace global;
void GeneratePartition(int partition_num) {
  int vertex_num = multi_tet->vertex_num_;
  std::ofstream out(dj::Format("%s/partition_tmp.txt", GetModelFolder()).c_str());
  out << partition_num << "\n";
  for (int i = 0; i < partition_num; ++i) {
    out << i << "\n";
  }
  out << vertex_num << "\n";
  int vnum_per_partition = (vertex_num - 1) / partition_num + 1;
  for (int i = 0; i < vertex_num; ++i) {
    // vertex_id partition_id local_vertex_id;
    out << i << " ";
    out << i / vnum_per_partition << " ";
    out << i % vnum_per_partition << "\n";
  }
  out.close();
}

void CreateBody() {
  L("Start creating bodies");
  std::vector<std::string>* transform_cmd = conf.Get<std::vector<std::string>* >("affine_transform");
  AffineTransformer<double> transform(*transform_cmd);
  //  body::body_mesh_file[0] = global::data_directory + body::body_mesh_file[0];
  //  auto mesh_file = dj::Format("%s/%s", GetModelFolder(), GetMeshFile());
  multi_tet = new MultiDomainTet(GetMeshFile(), 0, &transform);

  // TODO uncomment
  //  multi_tet->AttachSkeleton(global::skeleton);
  L("finish creating body");
  //  GeneratePartition(pose_conf.Get<int>("partition_num"));
  //    L("done"); exit(0);
}

void CreateSkeleton() {
  return;
  L("Create Skeleton");
  skeleton = new Skeleton(global::skeleton_file.c_str(), true);
}


void Init() {
  CreateSkeleton();
  CreateBody();
  global::current_body = multi_tet;
  //  multi_tet->LoadMeshPartition(DATA_DIRECTORY "armadillo/armadillo_150k.partition.txt");
  multi_tet->LoadPartitionInfo(GetPartitionFolder());
  //    multi_tet->LoadPartitionInfo(DATA_DIRECTORY "armadillo/one_partition.txt");
  //  multi_tet->LoadPartitionInfo(DATA_DIRECTORY "armadillo/partition_tmp.txt");
  //  multi_tet->LoadMass(DATA_DIRECTORY "armadillo/vertex_mass.txt");
  //  multi_tet->LoadLocalMass(DATA_DIRECTORY "armadillo//vertexMass");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/Basis2/Local_Modal_Basis");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/Local_Modal_Basis");
  //  multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/Basis2/Local_Pose_Basis");

  // FIXME uncomment
  // TODO uncomment
  //  L("No subspace loaded")
  auto GetSubspaceFileName = [](int part_id) -> const char* {
    static std::string name;
    name = dj::Format("%s/modal_basis/partition_%d.basis.bin", GetDataFolder(), part_id);
    return name.c_str();
  };
  multi_tet->LoadSubspace(GetSubspaceFileName, MultiDomainTet::kBinary);
  //  multi_tet->AssembleGlobalBasis();
  multi_tet->LoadCubature(dj::Format("%s/cubature", GetDataFolder()).c_str());
  multi_tet->PrecomputeFastSandwichTransform();

  std::vector<int> fixed_domain = *conf.Get<std::vector<int>*>("fixed domain");
  if (fixed_domain[0] != -1) {
    multi_tet->SetFixedDomains(std::set<int>(fixed_domain.begin(), fixed_domain.end()));
  }
  //  ASSERT(fixed_domain < multi_tet->part_num_);
  //  if (fixed_domain >= 0) {
  //    multi_tet->set_fixed_domain(fixed_domain);
  //  }

  int full_rigid_simulation = conf.Get<int>("full rigid motion");
  if (full_rigid_simulation) {
    multi_tet->EnableFullRigidSimulation();
  }

  if (conf.Get<int>("embed obj mesh")) {
    std::string obj_file = dj::Format("%z/../../3D Model of Octopus/OCTOPUS/OCTOPUS/octopus.3.obj", DATA_DIRECTORY);
    static ObjMesh obj_(obj_file.c_str(), ObjMesh::ASCII, false);
    static ObjMeshRender obj_render(&obj_);
    int textureMode = OBJMESHRENDER_GL_USEANISOTROPICFILTERING | OBJMESHRENDER_GL_USEMIPMAP | OBJMESHRENDER_LIGHTINGMODULATIONBIT;
    obj_render.loadTextures(textureMode);
//    multi_tet->UpdateObjMeshPosition();
    std::string mapping_folder = dj::Format("%z/obj_mapping", GetModelFolder());
    if (1) {
      multi_tet->EmbededObjMesh(&obj_, &obj_render, true, mapping_folder.c_str());
      multi_tet->EnalbeTexturedSurfaceRendering(DATA_DIRECTORY "../shader/octopus");
    } else {
      // resize mesh
      std::vector<double> verts;
      std::vector<int> tris;
      verts.resize(obj_.getNumVertices() * 3);
      tris.resize(obj_.getNumFaces() * 3);
      obj_.initTriangleLookup();
      for (int v = 0; v < int(obj_.getNumVertices()); ++v) {
        Vec3d pos = obj_.getPosition(v);
        verts[v * 3 + 0] = pos[0];
        verts[v * 3 + 1] = pos[1];
        verts[v * 3 + 2] = pos[2];
      }
      for (int t = 0; t < int(obj_.getNumFaces()); ++t) {
        obj_.getTriangle(t, &tris[t * 3 + 0], &tris[t * 3 + 1], &tris[t * 3 + 2]);
      }
      std::vector<std::string> cmd;
      cmd.push_back("centerize");
      cmd.push_back("range 1");
      cmd.push_back("min_y 0");
      cmd.push_back("rotate_y 180");
      AffineTransformer<double> transform(cmd);
      transform.Transform(&verts[0], int(verts.size() / 3));
      for (int v = 0; v < int(obj_.getNumVertices()); ++v) {
        Vec3d pos(verts[v * 3 + 0], verts[v * 3 + 1], verts[v * 3 + 2]);
        obj_.setPosition(v, pos);
        //        double dist = dj::Distance3(&pos[0], multi_tet->X + v * 3);
      }
      //      ObjMeshIO::Instance()->Write("/Users/dj/octopus.4.obj", verts, tris);
      //      exit(0);
      multi_tet->EmbededObjMesh(&obj_, &obj_render, false, mapping_folder.c_str());
      exit(0);
    }
    //  multi_tet->GenerateCubicaStylePartition(DATA_DIRECTORY "armadillo/cubica_partition"); L("cubica partition generated"); exit(0);
    //  multi_tet->SetConstrainedVertex();
    //    multi_tet->LoadSubspace(DATA_DIRECTORY "armadillo/basis/ascMat");
    //    multi_tet->SaveConstrainedVertex(DATA_DIRECTORY "armadillo/armadillo_150k.fixed_vertex.txt");
  } else {
    multi_tet->EnalbeSurfaceRenderingWithVBO(DATA_DIRECTORY "../shader/phong");
  }
  if (conf.Get<int>("export video")) {
    global::gl->EnableVideoRecording();
  }
  //  ASSERT(glGetError() == GL_NO_ERROR);
}
} // namespace pose_simulator

using namespace multi_domain_simulator;

int MultiDomainSimulator::HandleKeyPress(QKeyEvent *e) {
  switch (e->key()) {
    case Qt::Key_R:
      render_mode = (render_mode + 1) % MultiDomainTet::kRenderModeNum;
      break;
    case Qt::Key_X: {
      char file_name[512];
      sprintf(file_name, DATA_DIRECTORY "body%d.obj", 0);
      L(dj::Format("Saved body %z", file_name));
      L("Save body " + std::string(file_name));
      global::current_body->Save_OBJ_File(file_name);
      break;
    }
    default:
      return -1;
  }
  return kPriority_;
}

int MultiDomainSimulator::Render() {
  //  multi_tet->Render(MultiDomainTet::kRenderMode[render_mode], multi_tet->X);
  profiler.Start("render");
  multi_tet->Render(MultiDomainTet::kRenderMode[render_mode], multi_tet->X);
  profiler.End("render");
  return kNotHandled;
}

int MultiDomainSimulator::Idle() {
  static double total_simulation_time = conf.Get<double>("total simulation time");
  static unsigned int steps = 0;
  for (int i = 0; i < global::simulation_step_per_idle; ++i) {
    if (total_simulation_time > 0 && steps * global::time_step > total_simulation_time) {
      global::simulate = false;
      exit(0);
      continue;
    }
    steps++;
    //            multi_tet->FullSimulation(global::time_step);
    //            multi_tet->MultiBodyFullSimulation(global::time_step);
    //    multi_tet->MultiBodyFullSimulationExplicit(global::time_step);
    //        multi_tet->SubspaceMultiBodySimulation(global::time_step);
    //            multi_tet->SubspaceMultiBodySimulationWithCubature(global::time_step);
    //            multi_tet->SubspaceMultiBodySimulationOneStep(global::time_step);
    //    multi_tet->SubspaceMultiBodySimulationMultiIteration(global::time_step, 100);
    multi_tet->SubspaceMultiBodySimulationOneStepWithCubature(global::time_step);
    //  multi_tet->Simulate(global::time_step);
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
  //  multi_tet->UpdateSurfaceMeshVBO();
  //  ASSERT(glGetError() == GL_NO_ERROR);
  //  multi_tet->UpdateTexturedSurfaceMeshVBO();
  //  ASSERT(glGetError() == GL_NO_ERROR);
  return kNotHandled;
}

MultiDomainSimulator::MultiDomainSimulator()
  : InputHandler(2) {
  Init();
}

MultiDomainSimulator::~MultiDomainSimulator() {
  delete multi_tet;
}
