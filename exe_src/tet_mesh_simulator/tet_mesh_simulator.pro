#TARGET = main
# multi domain subspace simulation

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
include($$EXE_SRC_DIR/tet_mesh_simulator/tet_mesh_simulator.pri)
