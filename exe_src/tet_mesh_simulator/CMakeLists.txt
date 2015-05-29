set(sources
  tet_mesh_simulator_main.cpp
)

add_executable(tet_mesh_simulator ${headers} ${sources})
target_link_libraries(tet_mesh_simulator  ${LIBS})
