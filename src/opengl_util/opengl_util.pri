INCLUDEPATH += $$PWD
HEADERS += \
  $$PWD/camera.h \
  $$PWD/enviroment_cube_map.h \
  $$PWD/opengl_header.h \
  $$PWD/opengl_helper.h \
  $$PWD/opengl_soft_shadow_render.h \
  $$PWD/opengl_video_recorder.h \
  $$PWD/scene.h \
  $$PWD/textured_triangle_mesh_renderer.h \
  $$PWD/triangle_mesh_render.h \
  $$PWD/rainbow_color.h \
  $$PWD/ppm_io.h \
    $$PWD/reflective_object_renderer.h \
    $$PWD/opengl_material.h

SOURCES += \
  $$PWD/camera.cpp \
  $$PWD/opengl_helper.cpp \
  $$PWD/opengl_soft_shadow_render.cpp \
  $$PWD/opengl_video_recorder.cpp \
  $$PWD/scene.cpp \
  $$PWD/textured_triangle_mesh_renderer.cpp \
  $$PWD/triangle_mesh_render.cpp \
  $$PWD/rainbow_color.cpp \
  $$PWD/enviroment_cube_map.cpp \
    $$PWD/ppm_io.cpp \
    $$PWD/reflective_object_renderer.cpp \
    $$PWD/opengl_material.cpp
