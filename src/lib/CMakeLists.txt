set(headers
  affine_transformer.h
  blas_wrapper.h
  buffer_manager.h
  config_file.h
  conjugate_gradient_solver.h
  cubic_solver.h
  fixed_matrix.h
  fixed_matrix_utility.h
  fixed_vector.h
  fixed_vector_utility.h
  fps_calculator.h
  high_performance_timer.h
  macro_constant.h
  opengl_helper.h
  print_macro.h
  profiler.h
  progress_bar.h
  quaternion.h
  rainbow_color.h
  random.h
  scene.h
  singleton.h
  svd.h
  timer.h
  utility_function.h
  vector_io.h
  vector_lib.h
)

set(sources
  affine_transformer.cpp
  config_file.cpp
  conjugate_gradient_solver.cpp
  opengl_helper.cpp
  progress_bar.cpp
  quaternion.cpp
  rainbow_color.cpp
  scene.cpp
)

add_library(lib ${headers} ${sources})
