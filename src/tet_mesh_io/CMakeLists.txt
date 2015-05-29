set(headers
  cubica_tet_mesh_io.h
  msh_tet_mesh_io.h
  tet_gen_mesh_io.h
  tetrahedral_mesh_io.h
  vega_tet_mesh_io.h
)

set(sources
  cubica_tet_mesh_io.cpp
  msh_tet_mesh_io.cpp
  tet_gen_mesh_io.cpp
  tetrahedral_mesh_io.cpp
  vega_tet_mesh_io.cpp
)

add_library(tet_mesh_io ${headers} ${sources})
