set(headers
  bvh.h
  collision_util.h
  tet_collider.h
  util.h
)

set(sources
  bvh.cpp
  tet_collider.cpp
)

add_library(bvh ${headers} ${sources})
