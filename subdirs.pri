TEMPLATE = subdirs
include($$PWD/path_info.pri)

SUBDIRS += library
library.file = $$EXE_SRC_DIR/library/library.pro

exes = \
  multi_domain_subspace_simulator \
#  beam_generator \
#  tet_mesh_simulator \
  subspace_simulator \
  mesh_viewer \
  mixed_multi_domain \
  partition_info_generator \
  cubature \
  isosurface2obj \
#  multi_domain_basis \
  single_domain_basis \
  mesh_partitioner \
#  sample_pose \
#  subspace_pose_sample \
#  subspace_pose_viewer \
#  msh2tetgen \
#  random_subspace_pose_viewer \
#  cubica \
#  pose_viewer \
#  isosurface_stuffing \

#  obj_render \
#  rod \
#  solver_test \
#  pbd_rod \

for(exe, exes) {
  SUBDIRS += $${exe}
  $${exe}.file = $$EXE_SRC_DIR/$${exe}/$${exe}.pro
  $${exe}.depends = library
}


