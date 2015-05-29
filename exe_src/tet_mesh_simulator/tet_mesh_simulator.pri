#TARGET = tet_mesh_simulator
INCLUDEPATH += $$PWD
#LIBS += -L$$OUT_PWD/.. -lcommon_library
include(../../path_info.pri)
include($$SRC_DIR/common.pri)
SOURCES += \
    $$PWD/tet_mesh_simulator_main.cpp
