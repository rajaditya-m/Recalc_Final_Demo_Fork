#include "procedural_beam.h"
#include <fstream>
#include "print_macro.h"
#include "string_formatter.h"

int main(int argc, char *argv[]) {
  (void) argc;
  (void) argv;
  double edge_len = 0.05;
  int points_per_row = 8;
  int points_per_col = 4;
  int num_of_cross_section = 20;
  ProceduralBeam gen((points_per_row - 1) * edge_len,
                     (points_per_col - 1) * edge_len,
                     (num_of_cross_section - 1) * edge_len,
                     points_per_row,
                     points_per_col,
                     num_of_cross_section);
  gen.generateTetsProcedurally(dj::Format("/tmp/beam_%zx%zx%z", points_per_row, points_per_col, num_of_cross_section).c_str());
  if (1) {
    std::ofstream out("/tmp/vert_partition.txt");
    ASSERT(out.is_open());
    int num_of_cross_section_per_domain = 10;
    int domain_num = num_of_cross_section / num_of_cross_section_per_domain;
    int v_num = points_per_col * points_per_row * num_of_cross_section;
    out << v_num << "\n";
    for (int p = 0; p < domain_num; ++p) {
      for (int i = 0; i < points_per_col * points_per_row * num_of_cross_section_per_domain; ++i) {
        out << p << "\n";
      }
    }
    out.close();
  }
}
