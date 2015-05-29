set(sources
  sample_pose.cpp
)

add_executable(sample_pose ${headers} ${sources})
target_link_libraries(sample_pose  ${LIBS})
