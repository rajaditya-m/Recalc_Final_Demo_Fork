#include <sstream>
#include <algorithm>
#include "affine_transformer.h"
#include "vector_lib.h"
#include "print_macro.h"


template <class T>
void AffineTransformer<T>::Translate(T *position_array, int num, T dx, T dy, T dz) {
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    position_array[v3 + 0] += dx;
    position_array[v3 + 1] += dy;
    position_array[v3 + 2] += dz;
  }
}

template <class T>
void AffineTransformer<T>::GetCenter(const T *position_array, int num, T *center) {
  center[0] = center[1] = center[2] = 0;
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    center[0] += position_array[v3 + 0];
    center[1] += position_array[v3 + 1];
    center[2] += position_array[v3 + 2];
  }
  center[0] /= num;
  center[1] /= num;
  center[2] /= num;
}

template <class T>
void AffineTransformer<T>::Centerize(T *position_array, int num, T *center) {
  T original_center[3];
  GetCenter(position_array, num, original_center);
  if (center != NULL) {
    T translation[3] = {
      center[0] - original_center[0],
      center[1] - original_center[1],
      center[2] - original_center[2]
    };
    Translate(position_array, num, translation[0], translation[1], translation[2]);
  } else {
    Translate(position_array, num, -original_center[0], -original_center[1], -original_center[2]);
  }
}

template <class T>
void AffineTransformer<T>::SetMaxX(T *position_array, int num, T max_x) {
  SetExtremePosition(position_array, num, 0, max_x, true);
}

template <class T>
void AffineTransformer<T>::SetMaxY(T *position_array, int num, T max_y) {
  SetExtremePosition(position_array, num, 1, max_y, true);
}

template <class T>
void AffineTransformer<T>::SetMaxZ(T *position_array, int num, T max_z) {
  SetExtremePosition(position_array, num, 2, max_z, true);
}

template <class T>
void AffineTransformer<T>::SetMinX(T *position_array, int num, T min_x) {
  SetExtremePosition(position_array, num, 0, min_x, false);
}

template <class T>
void AffineTransformer<T>::SetMinY(T *position_array, int num, T min_y) {
  SetExtremePosition(position_array, num, 1, min_y, false);
}

template <class T>
void AffineTransformer<T>::SetMinZ(T *position_array, int num, T min_z) {
  SetExtremePosition(position_array, num, 2, min_z, false);
}

template <class T>
void AffineTransformer<T>::SetExtremePosition(T *position_array, int num, int xyz, T extreme_position, bool is_max) {
  T old_extreme_val;
  if (is_max) {
    old_extreme_val = -1e10;
    for (int i = 0; i < num; ++i) {
      old_extreme_val = std::max(position_array[i * 3 + xyz], old_extreme_val);
    }
  } else {
    old_extreme_val = 1e10;
    for (int i = 0; i < num; ++i) {
      old_extreme_val = std::min(position_array[i * 3 + xyz], old_extreme_val);
    }
  }
  T translation = extreme_position - old_extreme_val;
  for (int i = 0; i < num; ++i) {
    position_array[i * 3 + xyz] += translation;
  }
}


template <class T>
void AffineTransformer<T>::Scale(T *position_array, int num, T sx, T sy, T sz) {
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    position_array[v3 + 0] *= sx;
    position_array[v3 + 1] *= sy;
    position_array[v3 + 2] *= sz;
  }
}


template <class T>
void AffineTransformer<T>::ScaleAroundCenter(T *position_array, int num, T sx, T sy, T sz) {
  T center[3];
  GetCenter(position_array, num, center);
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    position_array[v3 + 0] = center[0] + (position_array[v3 + 0] - center[0]) * sx;
    position_array[v3 + 1] = center[1] + (position_array[v3 + 1] - center[1]) * sy;
    position_array[v3 + 2] = center[2] + (position_array[v3 + 2] - center[2]) * sz;
  }
}

template <class T>
void AffineTransformer<T>::SetMaxRange(T *position_array, int num, T max_range) {
  int xyz;
  T range = GetMaxRange(position_array, num, &xyz);
  T scale = max_range / range;
  Scale(position_array, num, scale, scale, scale);
}

template <class T>
void AffineTransformer<T>::GetRange(T *position_array, int num, T &min_x, T &min_y, T &min_z, T &max_x, T &max_y, T &max_z) {
  min_x = min_y = min_z = 1e10;
  max_x = max_y = max_z = -1e10;
  for (int i = 0, v3 = 0; i < num; ++i, v3 += 3) {
    min_x = std::min(min_x, position_array[v3 + 0]);
    max_x = std::max(max_x, position_array[v3 + 0]);
    min_y = std::min(min_y, position_array[v3 + 1]);
    max_y = std::max(max_y, position_array[v3 + 1]);
    min_z = std::min(min_z, position_array[v3 + 2]);
    max_z = std::max(max_z, position_array[v3 + 2]);
  }
}


template <class T>
void AffineTransformer<T>::GetRange(T *position_array, int num, T *min, T *max) {
  GetRange(position_array, num,
           min[0], min[1], min[2],
           max[0], max[1], max[2]);
}


template <class T>
T AffineTransformer<T>::GetMaxRange(T *position_array, int num, int* dim) {
  int tmp;
  if (dim == NULL) dim = &tmp;
  T min[3], max[3];
  GetRange(position_array, num, min[0], min[1], min[2], max[0], max[1], max[2]);
  T range[3] = {
    max[0] - min[0],
    max[1] - min[1],
    max[2] - min[2],
  };
  T max_range = range[0];
  (*dim) = 0;
  if (max_range < range[1]) {
    max_range = range[1];
    (*dim) = 1;
  }
  if (max_range < range[2]) {
    max_range = range[2];
    (*dim) = 2;
  }
  return max_range;
}


template <class T>
void AffineTransformer<T>::RotateAroundX(T *position_array, int num, T angle_in_degree) {
  T angle_in_radian = dj::Degree2Radian(angle_in_degree);
  for (int i = 0, v3 = 0; i < num; i++, v3 += 3) {
    T ty = position_array[v3 + 1] * cos(angle_in_radian) - position_array[v3 + 2] * sin(angle_in_radian);
    T tz = position_array[v3 + 1] * sin(angle_in_radian) + position_array[v3 + 2] * cos(angle_in_radian);
    position_array[v3 + 1] = ty;
    position_array[v3 + 2] = tz;
  }
}


template <class T>
void AffineTransformer<T>::RotateAroundY(T *position_array, int num, T angle_in_degree) {
  T angle_in_radian = dj::Degree2Radian(angle_in_degree);
  for (int i = 0, v3 = 0; i < num; i++, v3 += 3) {
    T tx = position_array[v3 + 0] *  cos(angle_in_radian) + position_array[v3 + 2] * sin(angle_in_radian);
    T tz = position_array[v3 + 0] * -sin(angle_in_radian) + position_array[v3 + 2] * cos(angle_in_radian);
    position_array[v3 + 0] = tx;
    position_array[v3 + 2] = tz;
  }
}


template <class T>
void AffineTransformer<T>::RotateAroundZ(T *position_array, int num, T angle_in_degree) {
  T angle_in_radian = dj::Degree2Radian(angle_in_degree);
  for (int i = 0, v3 = 0; i < num; i++, v3 += 3) {
    T tx = position_array[v3 + 0] * cos(angle_in_radian) - position_array[v3 + 1] * sin(angle_in_radian);
    T ty = position_array[v3 + 0] * sin(angle_in_radian) + position_array[v3 + 1] * cos(angle_in_radian);
    position_array[v3 + 0] = tx;
    position_array[v3 + 1] = ty;
  }
}

template <class T>
void AffineTransformer<T>::Rotate(T *position_array, int num, T *rotation_axis, T angle_in_degree) {
  T normalized_axis[3] = {rotation_axis[0], rotation_axis[1], rotation_axis[2]};
  dj::Normalize3(normalized_axis);
  T rotation_matrix[3][3];
  T angle_in_radian = dj::Degree2Radian(angle_in_degree);
  T cosine = cos(angle_in_radian);
  T sine = sin(angle_in_radian);
  rotation_matrix[0][0] = cosine + normalized_axis[0] * normalized_axis[0] * (1 - cosine);
  rotation_matrix[0][1] = normalized_axis[0] * normalized_axis[1] * (1 - cosine) - normalized_axis[2] * sine;
  rotation_matrix[0][2] = normalized_axis[0] * normalized_axis[2] * (1 - cosine) + normalized_axis[1] * sine;

  rotation_matrix[1][0] = normalized_axis[0] * normalized_axis[1] * (1 - cosine) + normalized_axis[2] * sine;
  rotation_matrix[1][1] = cosine + normalized_axis[1] * normalized_axis[1] * (1 - cosine);
  rotation_matrix[1][2] = normalized_axis[1] * normalized_axis[2] * (1 - cosine) - normalized_axis[0] * sine;

  rotation_matrix[2][0] = normalized_axis[0] * normalized_axis[2] * (1 - cosine) - normalized_axis[1] * sine;
  rotation_matrix[2][1] = normalized_axis[1] * normalized_axis[2] * (1 - cosine) + normalized_axis[0] * sine;
  rotation_matrix[2][2] = cosine + normalized_axis[2] * normalized_axis[2] * (1 - cosine);
  T tmp[3];
  for (int i = 0, v3 = 0; i < num; i++, v3 += 3) {
    tmp[0] = position_array[v3 + 0];
    tmp[1] = position_array[v3 + 1];
    tmp[2] = position_array[v3 + 2];
    dj::MulMatrix3x3Vec<T>((T (*)[3]) rotation_matrix, tmp, &position_array[v3]);
  }
}

template <class T>
void AffineTransformer<T>::AffineTransform(T *position_array, int num, const std::vector<std::string> &instruction) {
  for (unsigned i = 0; i < instruction.size(); ++i) {
    std::stringstream str(instruction[i]);
    std::string cmd;
    str >> cmd;
    if (cmd == "translate") {
      T dx, dy, dz;
      str >> dx >> dy >> dz;
      Translate(position_array, num, dx, dy, dz);
    } else if (cmd == "centerize") {
      T center[3];
      str >> center[0];
      if (str.good()) {
        str >> center[1] >> center[2];
        Centerize(position_array, num, center);
      } else {
        Centerize(position_array, num, NULL);
      }
    } else if (cmd == "max_x") {
      T max;
      str >> max;
      SetMaxX(position_array, num, max);
    } else if (cmd == "max_y") {
      T max;
      str >> max;
      SetMaxY(position_array, num, max);
    } else if (cmd == "max_z") {
      T max;
      str >> max;
      SetMaxZ(position_array, num, max);
    } else if (cmd == "min_x") {
      T min;
      str >> min;
      SetMinX(position_array, num, min);
    } else if (cmd == "min_y") {
      T min;
      str >> min;
      SetMinY(position_array, num, min);
    } else if (cmd == "min_z") {
      T min;
      str >> min;
      SetMinZ(position_array, num, min);
    } else if (cmd == "scale") {
      T scale[3];
      str >> scale[0] >> scale[1] >> scale[2];
      Scale(position_array, num, scale[0], scale[1], scale[2]);
    } else if (cmd == "scale_around_center") {
      T scale[3];
      str >> scale[0] >> scale[1] >> scale[2];
      ScaleAroundCenter(position_array, num, scale[0], scale[1], scale[2]);
    } else if (cmd == "rotate_x") {
      T angle;
      str >> angle;
      RotateAroundX(position_array, num, angle);
    } else if (cmd == "rotate_y") {
      T angle;
      str >> angle;
      RotateAroundY(position_array, num, angle);
    } else if (cmd == "rotate_z") {
      T angle;
      str >> angle;
      RotateAroundZ(position_array, num, angle);
    } else if (cmd == "rotate") {
      T axis[3];
      T angle;
      str >> axis[0] >> axis[1] >> axis[2] >> angle;
      Rotate(position_array, num, axis, angle);
    } else if (cmd == "range") {
      T range;
      str >> range;
      SetMaxRange(position_array, num, range);
    } else if (cmd == "none") {
    } else {
      std::cerr << CURRENT_LINE << " => invalid affine transformation command " << instruction[i] << std::endl;
    }
  }
}

template class AffineTransformer<float>;
template class AffineTransformer<double>;

void AffineTransformTest() {
  std::vector<std::string> cmd;
  cmd.push_back("translate 1 1 1");
  cmd.push_back("translate -1 -1 0");
  cmd.push_back("translate -0 -0 -1");
  cmd.push_back("scale 2 2 2");
  cmd.push_back("scale 1 0.5 1");
  cmd.push_back("scale 0.5 1.0 0.5");
  cmd.push_back("translate 0.5 0.5 0.5");
  cmd.push_back("scale_around_center 0.5 0.5 0.5");
  cmd.push_back("max_x 2");
  cmd.push_back("max_y 2");
  cmd.push_back("max_z 2");
  cmd.push_back("min_x 0.0");
  cmd.push_back("min_y 0.0");
  cmd.push_back("min_z 0.0");
  cmd.push_back("rotate_x 90");
  cmd.push_back("rotate_y 90");
  cmd.push_back("rotate_z 90");
  cmd.push_back("rotate 1 1 1 90");
}
