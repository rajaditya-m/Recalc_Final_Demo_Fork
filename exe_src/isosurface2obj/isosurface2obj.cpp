#include <Eigen/Dense>
#include <string>
#include <fstream>
#include "print_macro.h"
#include "string_formatter.h"
#include "config_file.h"
#include "global.h"
#include "CIsoSurface.h"
#include "hammock_level_set.h"

int main(void) {
  int nx = 5;
  int ny = 7;
  double radius = 0.0285;
  double height = 2.0;
  double width = (height * nx) / ny;
  width = 0.6;
  HammockLevelSet hammock(width, height, radius, nx, ny);

  double cell_length = 0.0183;
  double x_min = -width / 2 - radius - width * 1.5;
  double y_min = -height / 2 - radius - height * 1.5;
  double z_min = -1.5 * radius;
  int cell_x_num = int(-2.0 * x_min / cell_length) + 2;
  int cell_y_num = int(-2.0 * y_min / cell_length) + 2;
  int cell_z_num = int(-2.0 * z_min / cell_length) + 2;
  int xy_slice = cell_x_num * cell_x_num;
  (void) xy_slice;
  int total = cell_x_num * cell_y_num * cell_z_num;
  std::vector<double> level_set(total);
  for (int z = 0, idx = 0; z < cell_z_num; ++z) {
    for (int y = 0; y < cell_y_num; ++y) {
      for (int x = 0; x < cell_x_num; ++x, ++idx) {
        double pos[3] = {
          x_min + x * cell_length,
          y_min + y * cell_length,
          z_min + z * cell_length,
        };
        level_set[idx] = hammock(&pos[0]);
      }
    }
  }
  CIsoSurface<double> isosurface;
  P(cell_x_num, cell_y_num, cell_z_num);
  isosurface.GenerateSurface(&level_set[0], 0, cell_x_num - 1, cell_y_num - 1, cell_z_num - 1, cell_length, cell_length, cell_length);
#ifdef __APPLE__
  isosurface.ExportObjMesh("/Users/dj/isosurface.obj");
#elif defined(_WIN32)
  isosurface.ExportObjMesh("C:/tmp/isosurface.obj");
#else
  isosurface.ExportObjMesh("/home/dj/isosurface.obj");
#endif
}

