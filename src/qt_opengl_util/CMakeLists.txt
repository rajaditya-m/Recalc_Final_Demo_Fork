set(headers
  input_handler.h
  qt_gl_3d_scene_navigator.h
  qt_object_selector.h
  selectable_object.h
)

set(sources
  input_handler.cpp
  qt_gl_3d_scene_navigator.cpp
  qt_object_selector.cpp
)

add_library(qt_opengl_util ${headers} ${sources})
