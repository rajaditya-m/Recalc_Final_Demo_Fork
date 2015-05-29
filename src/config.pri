#TARGET = skeleton
TEMPLATE = app
QT += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
CONFIG -= app_bundle
CONFIG += opengl
CONFIG += console
#CONFIG -= flat
CONFIG(debug, debug | release) {
  MODE = "debug"
} else {
  MODE = "release"
  DEFINES += EIGEN_NO_DEBUG
}

include(../path_info.pri)
TMP_FILE_PREFIX = "$$OUT_PWD/../$$MODE"

win32|win64 {
  #TMP_FILE_PREFIX = "D:/tmp/$$PWD"
#  TMP_FILE_PREFIX = "C:/tmp/$$MODE"
  DEFINES += 'CURRENT_DIRECTORY=$$ROOT_DIR/data/'
  QMAKE_CXXFLAGS += /O2
} else {
#  TMP_FILE_PREFIX = "/tmp/$$ROOT_DIR/$$MODE"
  QMAKE_CXXFLAGS += -fopenmp
  QMAKE_CXXFLAGS += -Wno-unused-local-typedefs
  QMAKE_CXXFLAGS_RELEASE += -O3
  QMAKE_LFLAGS_RELEASE += -O3
  QMAKE_CXXFLAGS_RELEASE -= -O2
  QMAKE_LFLAGS_RELEASE -= -O1
  QMAKE_CXXFLAGS_DEBUG += -Og
  DEFINES += 'CURRENT_DIRECTORY=\'\"$$ROOT_DIR/data/\"\''
  DEFINES += 'DATA_DIRECTORY=\'\"$$ROOT_DIR/data/\"\''
}

LOCAL_LIB = "$${OUT_PWD}/../$${QMAKE_PREFIX_STATICLIB}common_library.$${QMAKE_EXTENSION_STATICLIB}"
PRE_TARGETDEPS += $${LOCAL_LIB}
LIBS += $${LOCAL_LIB}

macx {
  CC_VERSION = 4.8
#  QMAKE_CC = "/opt/local/bin/ccache /opt/local/bin/gcc-mp-$$CC_VERSION"
#  QMAKE_CXX = "/opt/local/bin/ccache /opt/local/bin/g++-mp-$$CC_VERSION"
  QMAKE_CC = "/opt/local/bin/gcc-mp-$$CC_VERSION"
  QMAKE_CXX = "/opt/local/bin/g++-mp-$$CC_VERSION"
  QMAKE_LINK = /opt/local/bin/g++-mp-$$CC_VERSION
  QMAKE_CXXFLAGS += -std=c++11
  QMAKE_CXXFLAGS += -Wno-deprecated
  QMAKE_CXXFLAGS += -Wno-deprecated-declarations
  QMAKE_LFLAGS += -fopenmp
  LIBS+= -framework GLUT -framework OpenGL -L/opt/local/lib -lglew
  LIBS += -lgslcblas
  INCLUDEPATH += /opt/local/include
  INCLUDEPATH += /usr/include
  INCLUDEPATH += /opt/local/include/eigen3
} else:unix {
#  ICC = /opt/intel/bin/icc
#  QMAKE_CC = $$ICC
#  QMAKE_CXX = $$ICC
#  QMAKE_LINK = $$ICC
  QMAKE_CC = "ccache gcc-4.8"
  QMAKE_CXX = "ccache g++-4.8"
  QMAKE_LINK = g++-4.8
  QMAKE_CXXFLAGS += -std=c++11  -Wno-unused-result
  INCLUDEPATH += /usr/include/eigen3
  LIBS += -lglut -lGL -lGLU -lGLEW -lgslcblas -lpthread -fopenmp
} else:win32|win64 {
  LIBS += -L"C:/lib/lib/x64"
  LIBS += arpack_x64.lib GlU32.lib glew32.lib freeglut_static.lib  OpenGL32.lib libgslcblas.a
  LIBS += -L"C:/Program Files (x86)/Intel/Composer XE 2015/compiler/lib/intel64"
  LIBS += -L"C:/Program Files (x86)/Intel/Composer XE 2015/mkl/lib/intel64"
  LIBS += mkl_intel_lp64.lib mkl_intel_thread.lib mkl_core.lib libiomp5md.lib
  INCLUDEPATH += "C:/Program Files (x86)/Intel/Composer XE 2015/mkl/include"
  INCLUDEPATH += "C:/lib/include"
}

msvc {
  #msvc:QMAKE_CXXFLAGS += -openmp
  QMAKE_CXXFLAGS += /openmp
  DEFINES += _CRT_SECURE_NO_WARNINGS
}

OTHER_FILES += \
    $$ROOT_DIR/conf/armadillo.conf \
    $$ROOT_DIR/conf/camera.conf \
    $$ROOT_DIR/shader/octopus.vert \
    $$ROOT_DIR/shader/octopus.frag \
    $$ROOT_DIR/shader/shadow.frag \
    $$ROOT_DIR/shader/shadow.vert \
    $$ROOT_DIR/shader/phong.vert \
    $$ROOT_DIR/shader/phong.frag \
    $$ROOT_DIR/shader/env_map.frag \
    $$ROOT_DIR/shader/env_map.vert \
    $$ROOT_DIR/CMakeLists.txt \

CONFIG += precompile_header

#message($${IN_PWD})
#message("+++++")
#message($$_PRO_FILE_)
#message($$_PRO_FILE_PWD_)
#message($${OUT_PWD})
#message($${HELLO})
#message("-------")
#DESTDIR = $$OUT_PWD
#OBJECTS_DIR = $$TMP_FILE_PREFIX/obj
#MOC_DIR = $$TMP_FILE_PREFIX/moc
#RCC_DIR = $$TMP_FILE_PREFIX/rcc
#UI_DIR = $$TMP_FILE_PREFIX/gui
#  DEFINES += USE_INTEL_MKL
#  INCLUDEPATH += /opt/intel/mkl/include
#  MKL_BLAS_LIB = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
#  MKL_BLAS_LIB_DIRECTORY = /opt/intel/mkl/lib/intel64
#  LIBS += -L$(MKL_BLAS_LIB_DIRECTORY) $(MKL_BLAS_LIB)



#message("$$IN_PWD")
