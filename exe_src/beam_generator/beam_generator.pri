INCLUDEPATH += $$PWD
include(../../path_info.pri)
include($$SRC_DIR/config.pri)
include($$SRC_DIR/common.pri)

SOURCES += \
  $$PWD/beam_generator.cpp

INCLUDEPATH += \
  $$SRC_DIR\lib \
