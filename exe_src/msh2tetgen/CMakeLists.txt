set(sources
  msh2tetgen.cpp
)

add_executable(msh2tetgen ${headers} ${sources})
target_link_libraries(msh2tetgen  ${LIBS})
