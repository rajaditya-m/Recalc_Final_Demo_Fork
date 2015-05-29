#TARGET = subspace_simulator
INCLUDEPATH += $$PWD
#LIBS += -L$$OUT_PWD/.. -lcommon_library
include(../../path_info.pri)
include($$SRC_DIR/common.pri)
#include($$SRC_DIR/basis_generation/basis_generation.pri)
SOURCES += \
    $$PWD/subspace_simulator_main.cpp
