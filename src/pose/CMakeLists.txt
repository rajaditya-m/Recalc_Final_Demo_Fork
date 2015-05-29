set(headers
  interactive_mesh_partitioner.h
  IO_FUNC.h
  MY_MATH.h
  pose_sampler.h
  pose_viewer_simulator.h
  sample_pose_simulator.h
  skeleton_adjuster.h
  skeleton_builder.h
  skeleton.h
  skeleton_node.h
  skeleton_rotator.h
  tet.h
  tet_mesh_simulator.h
  user_interaction.h
  VECTOR_TOOLS.h
)

set(sources
  interactive_mesh_partitioner.cpp
  pose_sampler.cpp
  pose_viewer_simulator.cpp
  sample_pose_simulator.cpp
  skeleton_builder.cpp
  skeleton.cpp
  skeleton_node.cpp
  skeleton_rotator.cpp
  tet.cpp
  tet_mesh_simulator.cpp
  user_interaction.cpp
)

add_library(pose ${headers} ${sources})
