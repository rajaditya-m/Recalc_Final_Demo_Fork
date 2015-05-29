set(headers
  multi_domain_subspace_simulator.h
  multi_domain_tet.h
  subspace_pose_sample_simulator.h
  subspace_pose_viewer.h
  subspace_simulator.h
  subspace_tet.h
)

set(sources
  multi_domain_subspace_simulator.cpp
  multi_domain_tet.cpp
  subspace_pose_sample_simulator.cpp
  subspace_pose_viewer.cpp
  subspace_simulator.cpp
  subspace_tet.cpp
)

add_library(subspace ${headers} ${sources})
