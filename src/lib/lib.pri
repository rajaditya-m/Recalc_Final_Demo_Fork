INCLUDEPATH += $$PWD

#HEADERS += \
PRECOMPILED_HEADER += \
    $$PWD/macro_constant.h \
    $$PWD/print_macro.h \
    $$PWD/fps_calculator.h \
    $$PWD/blas_wrapper.h \
    $$PWD/profiler.h \
    $$PWD/random.h \
    $$PWD/buffer_manager.h \
    $$PWD/timer.h \
    $$PWD/config_file.h \
    $$PWD/singleton.h \
    $$PWD/vector_lib.h \
    $$PWD/fixed_vector.h \
    $$PWD/utility_function.h \
    $$PWD/fixed_vector_utility.h \
    $$PWD/fixed_matrix_utility.h \
    $$PWD/fixed_matrix.h \
    $$PWD/affine_transformer.h \
    $$PWD/string_formatter.h \
    $$PWD/MERSENNETWISTER.h \
    $$PWD/matlab_io.h \
    $$PWD/conjugate_gradient_solver.h \
    $$PWD/quaternion.h \

HEADERS  += \
    $$PWD/triangle_mesh_creator.h \
    $$PWD/vector_io.h \
    $$PWD/svd.h \
    $$PWD/cubic_solver.h \
    $$PWD/progress_bar.h \
    $$PWD/biconjugate_gradient_stablized.h \
    $$PWD/statistics.h \
    $$PWD/text_file_io.h \
    $$PWD/binary_file_io.h \

SOURCES += \
    $$PWD/triangle_mesh_creator.cpp \
    $$PWD/conjugate_gradient_solver.cpp \
    $$PWD/config_file.cpp \
    $$PWD/affine_transformer.cpp \
    $$PWD/quaternion.cpp \
    $$PWD/progress_bar.cpp \
    $$PWD/biconjugate_gradient_stablized.cpp \
    $$PWD/string_formatter.cpp \
    $$PWD/text_file_io.cpp \
    $$PWD/binary_file_io.cpp \
    $$PWD/profiler.cpp \
   
