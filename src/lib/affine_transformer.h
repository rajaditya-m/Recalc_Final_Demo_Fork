#ifndef AFFINE_TRANSFORMER_H
#define AFFINE_TRANSFORMER_H
#include <vector>
#include <string>
template <class T = float>
class AffineTransformer {
public:
  static void Translate(T* position_array, int num, T dx, T dy, T dz);
  static void GetCenter(const T* position_array, int num, T* center);
  /// When center is NULL centerize at (0, 0, 0)
  static void Centerize(T *position_array, int num, T* center = NULL);

  static void SetMaxX(T* position_array, int num, T max_x);
  static void SetMaxY(T* position_array, int num, T max_y);
  static void SetMaxZ(T* position_array, int num, T max_z);
  static void SetMinX(T* position_array, int num, T min_x);
  static void SetMinY(T* position_array, int num, T min_y);
  static void SetMinZ(T* position_array, int num, T min_z);
  static void SetExtremePosition(T* position_array, int num, int xyz, T extreme_position, bool is_max);

  static void Scale(T* position_array, int num, T sx, T sy, T sz);
  static void ScaleAroundCenter(T* position_array, int num, T sx, T sy, T sz);

  static void SetMaxRange(T* position_array, int num, T max_range);
  static void GetRange(T* position_array, int num, T& min_x, T& min_y, T& min_z, T& max_x, T& max_y, T& max_z);
  static void GetRange(T* position_array, int num, T* min, T* max);
  static T GetMaxRange(T* position_array, int num, int *dim = NULL);

  static void RotateAroundX(T* position_array, int num, T angle_in_degree);
  static void RotateAroundY(T* position_array, int num, T angle_in_degree);
  static void RotateAroundZ(T* position_array, int num, T angle_in_degree);
  static void Rotate(T* position_array, int num, T* rotation_axis, T angle_in_degree);

  /// Available instructions:
  /// none
  /// translate dx dy dz
  /// centerize
  /// centerize x y z
  /// max_x max
  /// max_y max
  /// max_z max
  /// min_x min
  /// min_y min
  /// min_z min
  /// range max
  /// scale sx sy sz
  /// scale_around_center sx sy sz
  /// rotate_x angle_in_degree
  /// rotate_y angle_in_degree
  /// rotate_z angle_in_degree
  /// rotate rotation_axis angle_in_degree
  static void AffineTransform(T* position_array, int num, const std::vector<std::string>& instruction);

  AffineTransformer(std::vector<std::string>& transformation_command)
    : command_(transformation_command) {
  }

  void Transform(T* position_array, int num) {
    AffineTransform(position_array, num, command_);
  }

  std::vector<std::string> command_;
};

extern template class AffineTransformer<float>;
extern template class AffineTransformer<double>;

#endif // AFFINE_TRANSFORMER_H
