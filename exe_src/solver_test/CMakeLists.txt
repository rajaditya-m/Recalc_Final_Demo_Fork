set(sources
  solver_test.cpp
)

add_executable(solver_test ${headers} ${sources})
target_link_libraries(solver_test  ${LIBS})
