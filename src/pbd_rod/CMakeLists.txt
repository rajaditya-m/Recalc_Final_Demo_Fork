set(headers
  elastic_rod_shape.h
  frame_attachment.h
  pbd_elastic_rod.h
  pbd_rod_simulator.h
)

set(sources
  elastic_rod_shape.cpp
  frame_attachment.cpp
  pbd_elastic_rod.cpp
  pbd_rod_simulator.cpp
)

add_library(pbd_rod ${headers} ${sources})
