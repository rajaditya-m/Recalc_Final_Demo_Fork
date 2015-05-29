set(sources
  mesh_viewer_main.cpp
)

add_executable(mesh_viewer ${headers} ${sources})
target_link_libraries(mesh_viewer  ${LIBS})
