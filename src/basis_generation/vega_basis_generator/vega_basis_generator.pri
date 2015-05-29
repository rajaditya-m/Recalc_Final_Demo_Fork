INCLUDEPATH += $$PWD

HEADERS += \
  $$PWD/ARPACKSolver.h \
  $$PWD/CGSolver.h \
  $$PWD/computeStiffnessMatrixNullspace.h \
  $$PWD/cubicMesh.h \
  $$PWD/insertRows.h \
  $$PWD/invMKSolver.h \
  $$PWD/invKSolver.h \
  $$PWD/lapack-headers.h \
#  $$PWD/linearSolver.h \
  $$PWD/matrixBLAS.h \
  $$PWD/matrixMultiplyMacros.h \
  $$PWD/matrixPCA.h \
#  $$PWD/PardisoSolver.h \
#  $$PWD/sparseSolverAvailability.h \
  $$PWD/sparseSolvers.h \
  $$PWD/SPOOLESSolver.h \
  $$PWD/SPOOLESSolverMT.h \
  $$PWD/StVKCubeABCD.h \
  $$PWD/StVKElementABCDLoader.h \
  $$PWD/StVKHessianTensor.h \
  $$PWD/StVKTetHighMemoryABCD.h \
  $$PWD/triple.h \
  $$PWD/vegalong.h \
  $$PWD/performanceCounter.h \
  $$PWD/matrixBLASOptimized.cpp \
  $$PWD/matrixBLASVanilla.cpp \
  $$PWD/cubicMeshIntegrals.cpp \
  $$PWD/invKSolver.cpp \

SOURCES += \
  $$PWD/ARPACKSolver.cpp \
  $$PWD/CGSolver.cpp \
  $$PWD/computeStiffnessMatrixNullspace.cpp \
  $$PWD/cubicMesh.cpp \
  $$PWD/insertRows.cpp \
  $$PWD/invMKSolver.cpp \
  $$PWD/invKSolver.cpp \
#  $$PWD/linearSolver.cpp \
  $$PWD/matrixBLAS.cpp \
  $$PWD/matrixBLASOptimized.cpp \
  $$PWD/matrixBLASVanilla.cpp \
  $$PWD/matrixPCA.cpp \
#  $$PWD/PardisoSolver.cpp \
  $$PWD/SPOOLESSolver.cpp \
  $$PWD/SPOOLESSolverMT.cpp \
  $$PWD/StVKCubeABCD.cpp \
  $$PWD/StVKElementABCDLoader.cpp \
  $$PWD/StVKHessianTensor.cpp \
  $$PWD/StVKTetHighMemoryABCD.cpp \
