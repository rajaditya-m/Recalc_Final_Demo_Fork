DEFINES += VOLUMETRIC_ROD

INCLUDEPATH += $$PWD

HEADERS  += \
    $$PWD/tetrahedral_mesh.h \
    $$PWD/mass_spring_volumetric_object.h \
    $$PWD/volumetric_rod_simulator.h \
    $$PWD/fem_volumetric_object.h \
    $$PWD/subspace_mass_spring_volumetric_object.h \
    $$PWD/rigid_body.h \
    $$PWD/simulation_tetrahedral_mesh.h \

SOURCES += \
    $$PWD/tetrahedral_mesh.cpp \
    $$PWD/mass_spring_volumetric_object.cpp \
    $$PWD/volumetric_rod_simulator.cpp \
    $$PWD/fem_volumetric_object.cpp \
    $$PWD/subspace_mass_spring_volumetric_object.cpp \
    $$PWD/rigid_body.cpp \
    $$PWD/simulation_tetrahedral_mesh.cpp \
