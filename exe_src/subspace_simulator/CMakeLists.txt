set(sources
  subspace_simulator_main.cpp
)

add_executable(subspace_simulator ${headers} ${sources})
target_link_libraries(subspace_simulator  ${LIBS})
