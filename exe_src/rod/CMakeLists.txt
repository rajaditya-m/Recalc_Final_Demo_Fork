set(sources
  rod.cpp
)

add_executable(rod ${headers} ${sources})
target_link_libraries(rod  ${LIBS})
