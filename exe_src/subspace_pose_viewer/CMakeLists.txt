set(sources
  subspace_pose_viewer_main.cpp
)

add_executable(subspace_pose_viewer ${headers} ${sources})
target_link_libraries(subspace_pose_viewer  ${LIBS})
