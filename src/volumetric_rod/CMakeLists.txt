set(headers
  fem_volumetric_object.h
  mass_spring_volumetric_object.h
  rigid_body.h
  simulation_tetrahedral_mesh.h
  subspace_mass_spring_volumetric_object.h
  tetrahedral_mesh.h
  volumetric_rod_simulator.h
)

set(sources
  fem_volumetric_object.cpp
  mass_spring_volumetric_object.cpp
  rigid_body.cpp
  simulation_tetrahedral_mesh.cpp
  subspace_mass_spring_volumetric_object.cpp
  tetrahedral_mesh.cpp
  volumetric_rod_simulator.cpp
)

add_library(volumetric_rod ${headers} ${sources})
