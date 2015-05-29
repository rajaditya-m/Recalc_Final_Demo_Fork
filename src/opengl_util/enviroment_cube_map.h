#ifndef ENVIROMENT_CUBE_MAP_H_
#define ENVIROMENT_CUBE_MAP_H_
#include <string>

template <class Real>
class EnvironmentCubeMap {
  typedef unsigned int uint;
public:
  enum {
    kRealSize = sizeof(Real)
  };
  EnvironmentCubeMap(const char* cube_map_file);
  void RenderEvironmentCube(void);
  virtual ~EnvironmentCubeMap();

protected:
  float cube_vertecies_[6][4][3];
  float cube_tex_coordinate_[6][4][2];
  uint cube_map_face_[6];
  uint cube_map_handle_;
  void CreateCube(void);
};
#endif
