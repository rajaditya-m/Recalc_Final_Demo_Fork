include($$PWD/config.pri)

modules = \
  global \
  bvh \
  qt_opengl_util \
  lib \
  pose \
  vega \
  solver \
  volumetric_rod \
  window \
  pbd_rod \
  subspace \
  tet_mesh_io \
  tri_mesh_io \
  lma_calculator \
  cubature_optimization \
  cubica \
  procedural_beam \
  vega_obj_mtl \
  opengl_util \
  basis_generation \
  marching_cube \
  demo \
#  shape_op \
# cross-section-fem \
# obj_mtl_renderer \

for(module, modules) {
#  include($$PWD/$${module}/$${module}.pri)
  INCLUDEPATH += $$PWD/$$module
}

