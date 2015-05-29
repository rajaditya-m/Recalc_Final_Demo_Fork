#include "subspace_mass_spring_volumetric_object.h"
#include "rigid_body.h"
#include "volumetric_rod_simulator.h"
#include "mass_spring_volumetric_object.h"
#include "fem_volumetric_object.h"
#include "affine_transformer.h"
#include "tet_gen_mesh_io.h"
#include "opengl_helper.h"
#include "rainbow_color.h"
#include "qt_object_selector.h"
#include "open_gl_qt.h"

#include "tet_mesh_simulator_bridge.h"
#include "msh_tet_mesh_io.h"
TetMeshSimulatorBridge* vega;
class VolumetricRodSimulator::SimulatorInternal {
public:
  typedef FemVolumetricObject ElasticObj;
  //        typedef MassSpringVolumetricObject ElasticObj;
  //  typedef SubspaceMassSpringVolumetricObject ElasticObj;
  //      typedef RigidBody ElasticObj;


  SimulatorInternal() {
    //    Reduced();
    NonReducedOneObj();
    selector_ = new QtObjectSelector<Real>(rod_, Qt::LeftButton, Qt::ControlModifier);
    global::gl->InstallHandler(selector_);
  }


  void NonReducedOneObj() {
    //                std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/P";
    //                std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/staypuft";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/man";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/mid";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/low";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/dragon_vega";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/cube";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/two_tet";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/one_tet";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/tire";
    //        std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam_2seg";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam_20x20_1";
    //    std::string mesh_file = "/Users/dj/Downloads/global";
    std::string mesh_file = DATA_DIRECTORY "armadillo/armadillo_150k";
    //    std::string mesh_file = "/Users/dj/Dropbox/skeleton/data/armadillo/tmp/armadillo_400k.msh";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam";
    //                std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/hollow";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/hollow_single";
    //        std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/house2";

    std::vector<std::string> command;
    //    command.push_back("centerize");
    //    command.push_back("rotate_x -35");
    //    command.push_back("rotate_z 90");
    //        command.push_back("range 1.0");
    //    command.push_back("range 2.0");
    //    command.push_back("centerize 0 0.0 0");
    //    command.push_back("rotate_x -45");
    //    command.push_back("centerize");
    //    command.push_back("min_y 0.30000");
    //    command.push_back("min_x 0.00000");
    render_mode_ = 0;
    AffineTransformer<Real> trans0(command);
    rod_ = new ElasticObj(TetGenMeshIO::Instance(), mesh_file.c_str(), &trans0);
    //    rod_ = new ElasticObj(MshTetMeshIO::Instance(), mesh_file.c_str(), &trans0);
    //    vega = new TetMeshSimulatorBridge((SimulationTetrahedralMesh*) rod_);
    //    P(rod_->vert_[2][0]);
    //    rod_->vert_[2][0] = -rod_->vert_[2][0];
    P(rod_->v_num_);
    P(rod_->tet_num_);
  }

  void Reduced() {
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/hollow_single";
    //        std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/hollow8";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam_2seg";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam_20x40";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam_20x20_1";
    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam_20x20";
    //            std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/cube";
    //        std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/beam_small";
    //        std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/hollow";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/hollow_3seg";
    std::vector<std::string> command;
    command.push_back("centerize");
    //    command.push_back("rotate_x 45");
    //    command.push_back("range 2.0");
    command.push_back("centerize 0 0.0 0");
    command.push_back("min_y 0.300");
    //        command.push_back("min_x 0.300");
    //    command.push_back("min_z 0.300");
    render_mode_ = 0;
    AffineTransformer<Real> trans0(command);
    SubspaceMassSpringVolumetricObject* subspace_spring_obj = NULL;
    ASSERT(typeid(rod_) == typeid(subspace_spring_obj));
    subspace_spring_obj = new SubspaceMassSpringVolumetricObject(TetGenMeshIO::Instance(), mesh_file.c_str(), &trans0);
    //    subspace_spring_obj->LoadSubspace(DATA_DIRECTORY"subspace/basis-beam.txt");
    //    subspace_spring_obj->LoadSubspace(DATA_DIRECTORY"subspace/basis.txt");
    rod_ = (ElasticObj*) subspace_spring_obj;
    P(rod_->v_num_);
    P(rod_->tet_num_);
  }

  void NonReducedTwoObj() {
    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/P";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/staypuft";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/man";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/cube";
    //        std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/two_tet";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/one_tet";
    //        std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/tire";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/hollow";
    //    std::string mesh_file = DATA_DIRECTORY"tetgen_mesh/house2";
    std::vector<std::string> command;
    command.push_back("centerize");
    //    transform.push_back("rotate_x -90");
    command.push_back("range 0.5");
    command.push_back("centerize 0 0.0 0");
    command.push_back("min_y 0.00");
    AffineTransformer<Real> trans0(command);
    command.pop_back();
    command.push_back("centerize");
    //    transform.push_back("rotate_x -90");
    command.push_back("range 0.5");
    command.push_back("rotate_y 180");
    command.push_back("rotate_x 90");
    command.push_back("centerize 0 0.0 0");
    command.push_back("min_y 0.15");
    command.push_back("min_x -1.00");
    //        command.push_back("centerize 0.2 -0.1 -0.1");
    //    command.push_back("centerize -0.2 0 -0.1");
    //    command.push_back("centerize 0 0 0");
    //    command.push_back("rotate_x 90");
    //    command.push_back("scale 0.9 0.9 0.9");
    //    command.push_back("min_y 0.57");
    AffineTransformer<Real> trans1(command);
    render_mode_ = 0;
    //    render_mode_ = TetrahedraMesh::kWireFrame;
    std::vector<const char*> files(2, mesh_file.c_str());
    std::vector<AffineTransformer<Real>*> transforms;
    transforms.push_back(&trans0);
    transforms.push_back(&trans1);
    //    rod_ = new ElasticObj(TetGenMeshIO::Instance(), files, transforms);
    P(rod_->v_num_);
    P(rod_->tet_num_);
    //    rod_ = new ElasticObj(TetGenMeshIO::Instance(), mesh_file.c_str(), &trans0);

  }

  inline void Render() {
    rod_->Render(TetrahedralMesh::kRenderMode[render_mode_]);
  }

  ~SimulatorInternal() {
    delete rod_;
  }
  QtObjectSelector<Real>* selector_;
  ElasticObj* rod_;
  int render_mode_;
};


VolumetricRodSimulator::VolumetricRodSimulator()
  : InputHandler(1) {
  simulator_ = new SimulatorInternal();
}

VolumetricRodSimulator::~VolumetricRodSimulator() {
  delete simulator_;
}

int VolumetricRodSimulator::Idle() {
  profiler.Start("subspace");
  for (int i = 0; i < global::simulation_step_per_idle; ++i) {
    simulator_->rod_->Simulate(global::time_step);
    //    vega->Simulate(global::time_step);
  }
  profiler.End("subspace");
  return -1;
}

int VolumetricRodSimulator::Render() {
  simulator_->Render();
  return -1;
}

int VolumetricRodSimulator::HandleKeyPress(QKeyEvent *e) {
  if (e->key() == Qt::Key_R) {
    simulator_->render_mode_ = (simulator_->render_mode_ + 1) % TetrahedralMesh::kRenderModeNum;
  } else if (e->key() == Qt::Key_X) {
      char file_name[512];
      sprintf(file_name, DATA_DIRECTORY "body%d.obj", 0);
      simulator_->rod_->SaveSurface2Obj(file_name);
      L("Save body " + std::string(file_name));
  }
  return -1;
}
