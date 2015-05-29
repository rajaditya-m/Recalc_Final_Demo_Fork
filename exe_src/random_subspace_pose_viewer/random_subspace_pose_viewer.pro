#TARGET = mesh_partitioner
INCLUDEPATH += $$PWD
include(../../path_info.pri)
include($$SRC_DIR/common.pri)

PRO_FILE_NAME_FULL_PATH = $$_FILE_
PRO_FILE_BASE_NAME = $$basename(PRO_FILE_NAME_FULL_PATH)
MAIN_CPP = $$replace(PRO_FILE_BASE_NAME, .pro, .cpp)

PRO_FILE_NAME_FULL_PATH = $$_FILE_
PRO_FILE_BASE_NAME = $$basename(PRO_FILE_NAME_FULL_PATH)
MAIN_CPP = $$replace(PRO_FILE_BASE_NAME, .pro, .cpp)
EXE_NAME = $$replace(PRO_FILE_BASE_NAME, .pro, "")

win32|win64 {
 DEFINES += 'WINDOW_TITLE=$$EXE_NAME'
} else {
 DEFINES += 'WINDOW_TITLE=\'\"$$EXE_NAME\"\''
}

LIBS += -L$$OUT_PWD/.. -lcommon_library
SOURCES += \
  $$PWD/$$MAIN_CPP


#message($$PRO_FILE_NAME)
#message($$basename(PRO_FILE_NAME))
#message($$MAIN_CPP)
