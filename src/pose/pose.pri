INCLUDEPATH += $$PWD

SOURCES += \
    $$PWD/tet.cpp \
    $$PWD/skeleton.cpp \
    $$PWD/skeleton_node.cpp \
    $$PWD/user_interaction.cpp \
#    $$PWD/simulator.cpp \
    #$$PWD/skeleton_builder.cpp \
    $$PWD/skeleton_rotator.cpp \
    $$PWD/skeleton_builder.cpp \
    $$PWD/interactive_mesh_partitioner.cpp \
    $$PWD/pose_sampler.cpp \
    $$PWD/sample_pose_simulator.cpp \
    $$PWD/pose_viewer_simulator.cpp \
    $$PWD/tet_mesh_simulator.cpp

HEADERS  += \
    $$PWD/IO_FUNC.h \
    $$PWD/MY_MATH.h \
    $$PWD/skeleton_node.h \
    $$PWD/VECTOR_TOOLS.h \
    $$PWD/tet.h \
    $$PWD/skeleton.h \
    $$PWD/user_interaction.h \
#    $$PWD/simulator.h \
    $$PWD/skeleton_adjuster.h \
    $$PWD/skeleton_rotator.h \
    $$PWD/interactive_mesh_partitioner.h \
    $$PWD/pose_sampler.h \
    $$PWD/sample_pose_simulator.h \
    $$PWD/pose_viewer_simulator.h \
    $$PWD/tet_mesh_simulator.h
    #$$PWD/skeleton_builder.h \
