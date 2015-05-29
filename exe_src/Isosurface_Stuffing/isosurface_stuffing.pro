
TEMPLATE = app
TARGET = "Isosurface Stuffing"
INCLUDEPATH += . Util
QMAKE_CXXFLAGS += -std=c++11
macx {
  LIBS += -framework GLUT -framework OpenGL
} else:unix {
  LIBS += -lglut -lGLU
} else {
}
include($$PWD/../../path_info.pri)
include($$PWD/../../src/marching_cube/marching_cube.pri)
INCLUDEPATH += $$SRC_DIR
macx {
  INCLUDEPATH += /opt/local/include/eigen3
} else:unix {
  INCLUDEPATH += /usr/include/eigen3
} else:win32|win64 {
  INCLUDEPATH += "C:/lib/include"
  LIBS += -L"C:/lib/lib/x64"
  LIBS += GlU32.lib glew32.lib freeglut_static.lib  OpenGL32.lib libgslcblas.a
  INCLUDEPATH += "C:/lib/include"
}

# Input
HEADERS += Common.h \
           Edge.h \
           Mesh.h \
           Octant.h \
           Octree.h \
           Tetrahedron.h \
           Vertex.h \
           Util/Ball.h \
           Util/BallAux.h \
           Util/BallMath.h \
           Util/CoordSystem.h \
           Util/FrameSaver.h \
           Util/GLutilities.h \
           Util/inverse.h \
           Util/mathdefs.h \
           Util/Matrix.h \
           Util/MersenneTwister.h \
           Util/myMath.h \
           Util/noise.h \
           Util/Picker.h \
           Util/Quaternion.h \
           Util/SparseMatrix.h \
           Util/svd.h \
           Util/Timer.h \
           Util/Util.h \
           Util/vector.h \
           Util/vectorObj.h \
    net_level_set.h \
SOURCES += Edge.cpp \
           isosurface_stuffing.cpp \
           Mesh.cpp \
           Octant.cpp \
           Octree.cpp \
           Tetrahedron.cpp \
           Vertex.cpp \
           Util/BallAux.cxx \
           Util/BallMath.cxx \
           Util/CoordSystem.cxx \
           Util/FrameSaver.cpp \
           Util/GLutilities.cxx \
           Util/inverse.cpp \
           Util/Matrix.cpp \
           Util/myMath.cxx \
           Util/noise.cpp \
           Util/Picker.cpp \
           Util/Quaternion.cxx \
           Util/SparseMatrix.cpp \
           Util/svd.cpp \
           Util/Timer.cpp \
           Util/vector.cxx \
           Util/vectorObj.cxx
