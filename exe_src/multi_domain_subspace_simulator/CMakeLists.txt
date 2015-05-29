set(sources
  multi_domain_subspace_simulator_main.cpp
)

add_executable(multi_domain_subspace_simulator ${headers} ${sources})
target_link_libraries(multi_domain_subspace_simulator  ${LIBS})
