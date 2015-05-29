
PRO_FILE_NAME_FULL_PATH = $$_FILE_
PRO_FILE_BASE_NAME = $$basename(PRO_FILE_NAME_FULL_PATH)
MAIN_CPP = $$replace(PRO_FILE_BASE_NAME, .pro, .cpp)
EXE_NAME = $$replace(PRO_FILE_BASE_NAME, .pro, "")

win32|win64 {
 DEFINES += 'WINDOW_TITLE=$$EXE_NAME'
} else {
 DEFINES += 'WINDOW_TITLE=\'\"$$EXE_NAME\"\''
}
include(../../path_info.pri)
include($$EXE_SRC_DIR/mixed_multi_domain/mixed_multi_domain.pri)
include($$SRC_DIR/basis_generation/basis_generation.pri)
