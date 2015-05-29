include($$PWD/cubature_base/cubature_base.pri)

INCLUDEPATH += $$PWD
INCLUDEPATH += $$PWD/cubature_base

HEADERS += \
  $$PWD/CubatureObjDomain.h \
  $$PWD/CubatureObjInterface.h \
  $$PWD/multi_domain_cubature.h \
    $$PWD/single_domain_cubature.h

SOURCES += \
  $$PWD/CubatureObjInterface.cpp \
  $$PWD/CubatureObjDomain.cpp \
  $$PWD/multi_domain_cubature.cpp \
    $$PWD/single_domain_cubature.cpp
