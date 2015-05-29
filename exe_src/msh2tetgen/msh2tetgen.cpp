#include "tet_gen_mesh_io.h"
#include "msh_tet_mesh_io.h"
#include "global.h"
#include "tet.h"

int main(int argc, char *argv[]) {
  (void) argc; (void) argv;
  std::string msh = DATA_DIRECTORY "octopus/octopus.msh";
  std::string tetgen = DATA_DIRECTORY "octopus/octopus";
  std::vector<double> verts;
  std::vector<int> tets;
  MshTetMeshIO::Instance()->Read(msh.c_str(), verts, tets);
  TetGenMeshIO::Instance()->Write(tetgen.c_str(), verts, tets);
  Tet tet(tetgen.c_str(), 0, NULL);
  tet.Save_OBJ_File(DATA_DIRECTORY "octopus/octopus.obj");
}

