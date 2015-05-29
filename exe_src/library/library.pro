INCLUDEPATH += $$PWD
#IS_LIB = "true"
include(../../path_info.pri)
include($$SRC_DIR/common.pri)
for(module, modules) {
#  message(library.pro: module \"$$module\" is added to library)
  include($${SRC_DIR}/$${module}/$${module}.pri)
}
#message(library.pro: finish adding library)

LIBS = ""
PRE_TARGETDEPS = ""
TEMPLATE = lib
CONFIG += staticlib
DESTDIR = $$OUT_PWD/..
TARGET = common_library

SOURCES += \
    test.cpp \

HEADERS += \

