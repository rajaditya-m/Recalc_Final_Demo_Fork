set(sources
  subspace_pose_sample_main.cpp
)

add_executable(subspace_pose_sample ${headers} ${sources})
target_link_libraries(subspace_pose_sample  ${LIBS})
