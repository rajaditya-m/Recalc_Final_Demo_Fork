int dummy_int;
#ifdef MUL_TET_CPP
#include <Eigen/Geometry>
#include <fstream>
#include <set>
#include <unordered_set>
#include <cstdio>
#include "tet_gen_mesh_io.h"
#include "string_formatter.h"
#include "matlab_io.h"
#include "multi_domain_tet.h"
#include "conjugate_gradient_solver.h"
#include "biconjugate_gradient_stablized.h"
#include "global.h"
#include "vector_io.h"
#include "print_macro.h"
#include "tet_mesh_simulator_bridge.h"
#include "vector_io.h"
#include "config_file.h"
#include "BLOCK_MATRIX_GRAPH.h"

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// embeded mesh

#endif

//#define MIX_TET_CPP
#ifdef MIX_TET_CPP
#include <queue>
#include <functional>
#include <limits.h>
#include <utility>
#include "rainbow_color.h"
#include "global.h"
#include "opengl_helper.h"
#include "conjugate_gradient_solver.h"
#include "config_file.h"
#include "vector_lib.h"
#include "tet_mesh_simulator_bridge.h"
#include "isotropicHyperelasticFEM.h"
#include "mixed_multi_domain_tet.h"
#include "biconjugate_gradient_stablized.h"

#include "BLOCK_MATRIX_GRAPH.h"
#include "matlab_io.h"
#include "mixed_sparse_matrix.h"

#endif
