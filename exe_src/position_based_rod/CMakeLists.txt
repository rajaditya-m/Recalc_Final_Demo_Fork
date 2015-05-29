set(sources
  pbd_rod.cpp
)

add_executable(position_based_rod ${headers} ${sources})
target_link_libraries(position_based_rod  ${LIBS})
