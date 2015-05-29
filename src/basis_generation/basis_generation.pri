ADD_BASIS_GENERATION_MODULE = "false"
macx {
  ADD_BASIS_GENERATION_MODULE = "true"
  equals(ADD_BASIS_GENERATION_MODULE, true) {
    INCLUDEPATH += /opt/intel/mkl/include
    COMPOSER_XE_FOLDER = $$system(ls /opt/intel/ | grep 'composer_xe_20[0-9][0-9].[0-9]*.[0-9]*')
    LIBS += -L/opt/intel/$$COMPOSER_XE_FOLDER/compiler/lib
    LIBS += -L/opt/intel/mkl/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
    LIBS += -larpack
  }
} else:unix {
  ADD_BASIS_GENERATION_MODULE = "true"
  equals(ADD_BASIS_GENERATION_MODULE, true) {
    INCLUDEPATH += /opt/intel/mkl/include
    COMPOSER_XE_FOLDER = $$system(ls /opt/intel/ | grep 'composer_xe_20[0-9][0-9].[0-9]*.[0-9]*')
    LIBS += -L/opt/intel/$$COMPOSER_XE_FOLDER/compiler/lib/intel64
    LIBS += -L/opt/intel/mkl/lib/intel64
    LIBS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
    LIBS += -larpack
  }
}else{

  ADD_BASIS_GENERATION_MODULE = "true"
  equals(ADD_BASIS_GENERATION_MODULE, true) {
    LIBS += -L"C:/lib/lib/x64" arpack_x64.lib
  }
}
include($$PWD/pardio_solver/pardio_solver.pri);

#equals(ADD_BASIS_GENERATION_MODULE, true) {
  DEFINES += HAS_BASIS_GENERATION_MODULE
message("enabled-------------------");
  include($$PWD/vega_basis_generator/vega_basis_generator.pri)
  INCLUDEPATH += $$PWD

  HEADERS += \
    $$PWD/cubica_basis_generator.h \
    $$PWD/basis_generator.h \
    $$PWD/multi_domain_basis_generator.h \
    $$PWD/single_domain_basis_generator.h \
    $$PWD/MatrixOps.h \
    $$PWD/RedSVD_2.hpp
    #$$PWD/redsvd.hpp \
    #$$PWD/redsvdIncr.hpp \
    #$$PWD/util.hpp

  SOURCES += \
    $$PWD/basis_generator.cpp \
    $$PWD/cubica_basis_generator.cpp \
    $$PWD/multi_domain_basis_generator.cpp \
    $$PWD/single_domain_basis_generator.cpp \
    $$PWD/MatrixOps.cpp
    #$$PWD/util.cpp
#}

