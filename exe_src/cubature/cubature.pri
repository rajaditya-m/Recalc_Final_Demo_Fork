#TARGET = subspace_simulator
INCLUDEPATH += $$PWD
#LIBS += -L$$OUT_PWD/.. -lcommon_library
include(../../path_info.pri)
include($$SRC_DIR/common.pri)
SOURCES += \
    $$PWD/cubature.cpp
