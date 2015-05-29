SET(exes
  mesh_viewer
  msh2tetgen
  multi_domain_subspace_simulator
  pose_viewer
  position_based_rod
  rod
  sample_pose
  solver_test
  subspace_pose_sample
  subspace_pose_viewer
  subspace_simulator
  tet_mesh_simulator
)


foreach(exe ${exes})
  add_subdirectory(${exe})
endforeach(exe)
	
