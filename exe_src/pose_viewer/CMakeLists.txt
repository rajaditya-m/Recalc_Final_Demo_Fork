set(sources
  pose_viewer.cpp
)

add_executable(pose_viewer ${headers} ${sources})
target_link_libraries(pose_viewer  ${LIBS})
