INCLUDEPATH += $$PWD
include(../../path_info.pri)
include($$SRC_DIR/solver/solver.pri)
#include($$SRC_DIR/lib/lib.pri)

SOURCES += \
  $$PWD/solver_test.cpp

INCLUDEPATH += \
  $$SRC_DIR\lib \
