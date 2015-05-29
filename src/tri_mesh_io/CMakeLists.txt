set(headers
  mesh_viewer.h
  obj_mesh_io.h
  off_mesh_io.h
  triangular_mesh.h
  triangular_mesh_io.h
)

set(sources
  mesh_viewer.cpp
  obj_mesh_io.cpp
  off_mesh_io.cpp
  triangular_mesh.cpp
  triangular_mesh_io.cpp
)

add_library(tri_mesh_io ${headers} ${sources})
target_link_libraries(tri_mesh_io tet_mesh_io)
