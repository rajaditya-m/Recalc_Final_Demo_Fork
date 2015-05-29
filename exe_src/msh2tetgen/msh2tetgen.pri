INCLUDEPATH += $$PWD
#LIBS += -L$$OUT_PWD/.. -lcommon_library
include(../../path_info.pri)
include($$SRC_DIR/common.pri)
#include($$SRC_DIR/tet_mesh_io/li.pri)

SOURCES += \
  $$PWD/msh2tetgen.cpp
