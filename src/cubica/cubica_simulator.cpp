#include "cubica_simulator.h"
#include "global.h"
#include "cubica_tet.h"
#include "open_gl_qt.h"
#include "qt_object_selector.h"
#include "skeleton_node.h"
#include "skeleton.h"
#ifdef HAS_BASIS_GENERATION_MODULE
#include "cubica_basis_generator.h"
#endif

namespace cubica_simulator {
CubicaTet* cubica_tet;
void Init() {
  // Generate basis
  if (0) {
#ifdef HAS_BASIS_GENERATION_MODULE
    CubicaBasisGenerator gen(DATA_DIRECTORY "armadillo/armadillo_150k",
                             DATA_DIRECTORY "armadillo/cubica_partition");
    global::skeleton = new Skeleton(global::skeleton_file.c_str(), true);
    gen.AttachSkeleton(global::skeleton);
    gen.GenerateBasis(DATA_DIRECTORY "armadillo/cubica_partition", 26, 20);
    L("finishing generating basis for cubica tet");
    exit(0);
#endif
  }
  cubica_tet = new CubicaTet(DATA_DIRECTORY "armadillo/armadillo_150k",
                             DATA_DIRECTORY "armadillo/cubica_partition");
  auto GetSubspaceFile = [](int part_id) -> const char* {
    static char file_name[512];
    sprintf(file_name, DATA_DIRECTORY "armadillo/cubica_partition/partition_%d.basis.bin", part_id);
    return file_name;
  };
  cubica_tet->LoadSubspace(GetSubspaceFile);
  cubica_tet->LoadCubature(GetSubspaceFile);
  global::gl->InstallHandler(new QtObjectSelector<double>(cubica_tet));
}

}

using namespace cubica_simulator;

CubicaSimulator::CubicaSimulator() {
  Init();
}

CubicaSimulator::~CubicaSimulator() {
  delete cubica_tet;
}

int CubicaSimulator::Render() {
  cubica_tet->Render();
  return kNotHandled;
}

int CubicaSimulator::Idle() {
  cubica_tet->Simulate(global::time_step);
  return kNotHandled;
}

int CubicaSimulator::HandleKeyPress(QKeyEvent *e) {
  switch (e->key()) {
    case Qt::Key_R:
      cubica_tet->NextRenderMode();
      break;
    default:
      return -1;
  }
  return kPriority_;
}
