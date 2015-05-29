#include "subspace_simulator.h"
#include "affine_transformer.h"
#include "tet_mesh_simulator_bridge.h"
#include "tet.h"
#include "input_handler.h"
#include "open_gl_qt.h"
#include "subspace_tet.h"
#include "string_formatter.h"
#include "vega_tet_mesh_io.h"
#include "qt_object_selector.h"
#include "global.h"
#include "skeleton.h"
#include "skeleton_node.h"
#include "pose_sampler.h"
#include "tet_gen_mesh_io.h"
#include "rainbow_color.h"
#include "opengl_header.h"
#include "config_file.h"

#include <fstream>

#include "basis_generator.h"


namespace subspace_simulator {
SubspaceTet* tet;
bool newBasisConstructed = false;
SingleDomainBasisGenerator *basisGenerator;
int part_id = 16;
int render_mode = 0;
using namespace global;
void CreateBody() {
  L("Start creating bodies");
  std::vector<std::string>* transform_cmd = conf.Get<std::vector<std::string>* >("affine_transform");
  AffineTransformer<double> transform(*transform_cmd);
  char file_name[512];
  sprintf(file_name, "%s", GetMeshFile(), part_id);
  tet = new SubspaceTet(file_name, 0, &transform);
  L("finish creating body");
}

void CreateSkeleton() {
  L("Create Skeleton");
  skeleton = new Skeleton(global::skeleton_file.c_str(), true);
}


void Init() {
  //CreateSkeleton();
  CreateBody();

  //Compute the basis first
  //SingleDomainBasisGenerator generator(GetMeshFile());
  std::vector<std::string>* transform_cmd = conf.Get<std::vector<std::string>* >("affine_transform");
  AffineTransformer<double> transform(*transform_cmd);
  basisGenerator = new SingleDomainBasisGenerator(GetMeshFile(),&transform);
  std::string basis_prefix = dj::Format("%z/modal_basis/genBasis",GetDataFolder());
  //std::string mass_file = dj::Format("%z/vertexmass.M",GetDataFolder());
  //basisGenerator->LoadMass(mass_file.c_str());
  std::string fixed_vertex_file = dj::Format("%z/fixed_verts.bou",GetModelFolder());
  basisGenerator->ProcessFixedVertex(fixed_vertex_file.c_str());
  basisGenerator->preLoad(basis_prefix.c_str());
  //basisGenerator->GenerateBasis(basis_prefix.c_str(),30,50);

  //Compute the cubature here
/*
  float scale = 0.01000;
  int maxCubaturePoints = 140;
  SingleDomainCubature singleDomainCubOp(GetMeshFile(),&transform);
  std::string output_folder = dj::Format("%z/cubature", GetDataFolder());
  std::string name = dj::Format("%s/modal_basis/genBasis.basis.bin", GetDataFolder());
  singleDomainCubOp.LoadBinarySubspace(name.c_str());
  singleDomainCubOp.SetFolder(100,output_folder,scale,true);//true);
  singleDomainCubOp.GenerateCubature(maxCubaturePoints, 1.0e-6);
*/

  /*float scale = 0.01000;
  int maxCubaturePoints = 90;
  std::string meshname = dj::Format("%s/stem", GetDataFolder());
  SingleDomainCubature singleDomainCubOp(meshname.c_str(),&transform,true,false,0);
  std::string output_folder = dj::Format("%z/cubature_1", GetDataFolder());
  std::string name = dj::Format("%s/modal_basis/genBasis.basis_1.bin", GetDataFolder());
  singleDomainCubOp.LoadBinarySubspace(name.c_str());
  singleDomainCubOp.SetFolder(100,output_folder,scale,true);//true);
  singleDomainCubOp.GenerateCubature(maxCubaturePoints, 1.0e-6);

  maxCubaturePoints = 170;
  meshname = dj::Format("%s/head", GetDataFolder());
  SingleDomainCubature singleDomainCubOp2(meshname.c_str(),&transform,true,false,0);
  output_folder = dj::Format("%z/cubature_2", GetDataFolder());
  name = dj::Format("%s/modal_basis/genBasis.basis_2.bin", GetDataFolder());
  singleDomainCubOp2.LoadBinarySubspace(name.c_str());
  singleDomainCubOp2.SetFolder(100,output_folder,scale,true);//true);
  singleDomainCubOp2.GenerateCubature(maxCubaturePoints, 1.0e-6);*/


  //This is the multidomain cubature merge rioutine
  std::vector<int> tid;
  std::vector<double> weights;
  int cubature_flip_counter_ = 0;
  int elem_offset = 29714;
  char file_name_2[512];
  sprintf(file_name_2, "%s/cubature_1/cubature_pnts.txt", GetDataFolder());
  std::ifstream in(file_name_2);
  ASSERT(in.is_open(), P(file_name_2));
  int cubature_point_num;
  in >> cubature_point_num;
  //  P(cubature_point_num);
  for (int i = 0; i < cubature_point_num; ++i) {
    int tet = -1;
    double weight = -1e20;
    in >> tet;
    in >> weight;
    if (weight < 1e-10) {
      continue;
      KK;
    } else {
      tid.push_back(tet);
      weights.push_back(weight);
      cubature_flip_counter_++;
    }
  }
  in.close();
  P(cubature_flip_counter_);
  tet->setCubatureFlipCounter(cubature_flip_counter_);
  sprintf(file_name_2, "%s/cubature_2/cubature_pnts.txt", GetDataFolder());
  in.open(file_name_2);
  ASSERT(in.is_open(), P(file_name_2));
  in >> cubature_point_num;
  //  P(cubature_point_num);
  for (int i = 0; i < cubature_point_num; ++i) {
    int tet = -1;
    double weight = -1e20;
    in >> tet;
    in >> weight;
    if (weight < 1e-10) {
      continue;
      KK;
    } else {
      tid.push_back(tet+elem_offset);
      weights.push_back(weight);
    }
  }
  in.close();
  sprintf(file_name_2, "%s/cubature/cubature_pnts.txt", GetDataFolder());
  std::ofstream cub_out(file_name_2);
  ASSERT(cub_out.is_open(), P(file_name_2));
  cub_out << tid.size() << std::endl;
  for (int i = 0; i < int(tid.size()); ++i) {
      cub_out << tid[i] << " " << weights[i] << std::endl;
  }
  cub_out.close();

  global::current_body = tet;
  char file_name[512];
  int numRows = tet->vertex_num_ * 3;
  int numBasis = basisGenerator->basis_generator->non_linear_mode_num_;
  tet->LoadBinarySubspace(basisGenerator->basis_generator->non_linear_modes_,numRows,numBasis);
  //int numBasis = basisGenerator->basis_generator->linear_mode_num_;
  //tet->LoadBinarySubspace(basisGenerator->basis_generator->pure_eigen_vectors_,numRows,numBasis);

  //sprintf(file_name, "%s/modal_basis/genBasis.basis.bin", GetDataFolder());
  //tet->LoadBinarySubspace(file_name);
  //tet->LoadCubature(singleDomainCubOp.getCubaturePoints(),singleDomainCubOp.getCubatureWeights());
  tet->LoadCubature(file_name_2);

  newBasisConstructed = false;

   tet->EnalbeSurfaceRenderingWithVBO(DATA_DIRECTORY "../shader/phong");

}

}

using namespace subspace_simulator;

SubspaceSimulator::SubspaceSimulator()
  : InputHandler(10) {
  Init();
  global::gl->InstallHandler(new QtObjectSelector<double>(tet));
}

SubspaceSimulator::~SubspaceSimulator() {
  delete tet;
}

int SubspaceSimulator::HandleKeyPress(QKeyEvent *e) {
  switch (e->key()) {
    case Qt::Key_R:
      render_mode = (render_mode + 1) % Tet::kRenderModeNum;
      break;
    case Qt::Key_X: {
      char file_name[512];
      sprintf(file_name, "%s/body%d.obj", GetModelFolder(), 0);
      L("Save body " + std::string(file_name));
      global::current_body->Save_OBJ_File(file_name);
      break;
    }
    case Qt::Key_K:
      //tet->NextBasis();
      global::ap1 += 0.0001;
      break;
    case Qt::Key_J:
      //tet->PrevBasis();
      global::ap1 -= 0.0001;
      if(global::ap1<0)
          global::ap1 = 0.001;
      break;
  case Qt::Key_V:
   global::sim_state = 1;
   tet->setSimMode(1);
   tet->Reset();
   //tet->saveOldBasis();
   tet->LoadBinarySubspace(basisGenerator->basis_generator->non_linear_modes_,tet->vertex_num_ * 3,basisGenerator->basis_generator->non_linear_mode_num_);
   newBasisConstructed = false;
   break;
  case Qt::Key_C:
    //tet->NextBasis();
    //global::ap2 += 0.005;
      global::stiffyMult *= 10;
      basisGenerator->setStiffyMult(global::stiffyMult);
    break;
  case Qt::Key_M:
    //tet->NextBasis();
    global::ap2 += 0.005;
      //global::stiffyMult *= 10;
      //basisGenerator->setStiffyMult(global::stiffyMult);
    break;
  case Qt::Key_N:
    //tet->PrevBasis();
    global::ap2 -= 0.00005;
    if(global::ap2<0)
        global::ap2 = 0.0005;
    break;
  case Qt::Key_S:
    //if(tet->getSimMode()==1)
    tet->setSimMode(2);
    //else
    //    tet->setSimMode(4);
    //tet->saveOldBasis();
    newBasisConstructed = false;
    break;
  case Qt::Key_G:
      P(global::ap1);
      P(global::ap2);
      break;
  default:
      return -1;
  }
  return kPriority_;
}

int SubspaceSimulator::Render() {
  tet->Render(Tet::kRenderMode[render_mode], tet->render_X);
  glColor3fv(kYellow());
  glBegin(GL_POINTS);
  Vertex3v(&(tet->center_of_mass_[0]));
  glEnd();
  return -1;
}

int SubspaceSimulator::Idle() {
    tet->SimulateWithCubature(global::time_step);
    //tet->SimulateWithRigidMotion(global::time_step);
    //tet->Simulate(global::time_step);
    int statusCheck = tet->getSimMode();
    if(statusCheck==2 && !newBasisConstructed) {
        global::sim_state = 2;
        global::gl->updateGL();
        newBasisConstructed = true;
        basisGenerator->setStitchedStiffnessMatrix(tet->inv_fem_->GetStitchStiffnessMatrixPointer(2));
          profiler.Start("Basis Gen");
        basisGenerator->RegenerateStitchBasis(2);
         profiler.End("Basis Gen");
         //tet->LoadBinarySubspace(basisGenerator->basis_generator->non_linear_modes_,tet->vertex_num_ * 3,basisGenerator->basis_generator->non_linear_mode_num_);
         tet->LoadBinarySubspace(basisGenerator->basis_generator->stitched_non_linear_modes_,(tet->vertex_num_*3),basisGenerator->basis_generator->stitched_non_linear_mode_num_);
         //tet->LoadBinarySubspace(basisGenerator->basis_generator->stitched_linear_modes_,(tet->vertex_num_*3),basisGenerator->basis_generator->stitched_linear_mode_num_);
        global::sim_state = 3;
    }
    else if(statusCheck==4 && !newBasisConstructed) {
        //P("Does it??");
        global::sim_state = 4;
        global::gl->updateGL();
        newBasisConstructed = true;
        basisGenerator->setStitchedStiffnessMatrix(tet->inv_fem_->GetStitchStiffnessMatrixPointer(4));
        basisGenerator->RegenerateStitchBasis(4);
        tet->LoadBinarySubspace(basisGenerator->basis_generator->stitched_non_linear_modes_,(tet->vertex_num_*3),basisGenerator->basis_generator->stitched_non_linear_mode_num_);
        global::sim_state = 5;
    }
    //tet->Save_OBJ_File();
  return -1;
}
