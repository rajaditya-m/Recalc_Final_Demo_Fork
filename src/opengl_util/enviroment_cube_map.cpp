#include <vector>
#include "enviroment_cube_map.h"
#include "ppm_io.h"
#include "opengl_header.h"


template <class Real>
EnvironmentCubeMap<Real>::EnvironmentCubeMap(const char *cube_map_file) {
  std::vector<unsigned char> cube_map_data[6];
  int width, height;
  std::string file_names[6];
  file_names[0] = std::string(cube_map_file) + "_posx.ppm";
  file_names[1] = std::string(cube_map_file) + "_negx.ppm";
  file_names[2] = std::string(cube_map_file) + "_posy.ppm";
  file_names[3] = std::string(cube_map_file) + "_negy.ppm";
  file_names[4] = std::string(cube_map_file) + "_posz.ppm";
  file_names[5] = std::string(cube_map_file) + "_negz.ppm";
  for (int i = 0; i < 6; ++i) {
    dj::ReadPPM(file_names[i].c_str(), cube_map_data[i], width, height) ;
  }
  glGenTextures(1, &cube_map_handle_);
  //    glActiveTexture(GL_TEXTURE0 + cubemap_texture);
  //    glEnable(GL_TEXTURE_CUBE_MAP);
  glBindTexture(GL_TEXTURE_CUBE_MAP, cube_map_handle_);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  for (int i = 0; i < 6; i++) {
    glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, &cube_map_data[i][0]);
  }

  for (int i = 0; i < 6; i++) {
    glGenTextures(1, &cube_map_face_[i]);
    glBindTexture(GL_TEXTURE_2D, cube_map_face_[i]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, width, height, 0,  GL_RGB, GL_UNSIGNED_BYTE, &cube_map_data[i][0]);
  }
  CreateCube();
}

template <class Real>
void EnvironmentCubeMap<Real>::RenderEvironmentCube() {
  glActiveTexture(GL_TEXTURE0);
  glDisable(GL_LIGHTING);
  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  for (int i = 0; i < 6; ++i) {
    glBindTexture(GL_TEXTURE_2D, cube_map_face_[i]);
    glBegin(GL_TRIANGLE_FAN);
    for (int j = 0; j < 4; j++) {
      glTexCoord2fv(&cube_tex_coordinate_[i][j][0]);
      glVertex3fv(&cube_vertecies_[i][j][0]);
    }
    glEnd();
  }
  glDisable(GL_TEXTURE_2D);
  //------------------------------------------------------------------------------
  //    glDisable(GL_LIGHTING);
  //    glEnable(GL_TEXTURE_2D);
  //    glBindTexture(GL_TEXTURE_2D, cube_map_face_[3]);
  //    glBegin(GL_QUADS);
  //    glTexCoord2f(0, 0);
  //    glVertex2f(0, 0);
  //    glTexCoord2f(0, 1);
  //    glVertex2f(0, 1);
  //    glTexCoord2f(1, 1);
  //    glVertex2f(1, 1);
  //    glTexCoord2f(1, 0);
  //    glVertex2f(1, 0);
  //    glEnd();

  //    glVertexPointer(3, GL_FLOAT, 0, &cube_vertecies_[0][0][0]);
  //    glTexCoordPointer(2, GL_FLOAT, 0, &tex_coordinate_[0][0][0]);
  //    glEnableClientState(GL_VERTEX_ARRAY);
  //    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
  //    glDisable(GL_LIGHTING);
  //    glEnable(GL_TEXTURE_2D);
  //    for (int i = 0; i < 6; ++i) {
  //      glBindTexture(GL_TEXTURE_2D, cube_map_face_[i]);
  //      glDrawArrays(GL_TRIANGLE_FAN, i * 4, 4);
  //    }
  //    glDisableClientState(GL_VERTEX_ARRAY);
  //    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
  //    glDisable(GL_TEXTURE_2D);
}

template <class Real>
EnvironmentCubeMap<Real>::~EnvironmentCubeMap() {
}


template <class Real>
void EnvironmentCubeMap<Real>::CreateCube() {
  //    static const float coords[6][4][3] = {
  //      { { +1, -1, -1 }, { -1, -1, -1 }, { -1, +1, -1 }, { +1, +1, -1 } }, // back
  //      { { +1, +1, -1 }, { -1, +1, -1 }, { -1, +1, +1 }, { +1, +1, +1 } },
  //      { { +1, -1, +1 }, { +1, -1, -1 }, { +1, +1, -1 }, { +1, +1, +1 } },
  //      { { -1, -1, -1 }, { -1, -1, +1 }, { -1, +1, +1 }, { -1, +1, -1 } },
  //      { { +1, -1, +1 }, { -1, -1, +1 }, { -1, -1, -1 }, { +1, -1, -1 } },
  //      { { -1, -1, +1 }, { +1, -1, +1 }, { +1, +1, +1 }, { -1, +1, +1 } }
  //    };
#define VERTEX_0 {-1, -1, -1}
#define VERTEX_1 {-1, -1, +1}
#define VERTEX_2 {+1, -1, +1}
#define VERTEX_3 {+1, -1, -1}
#define VERTEX_4 {-1, +1, -1}
#define VERTEX_5 {-1, +1, +1}
#define VERTEX_6 {+1, +1, +1}
#define VERTEX_7 {+1, +1, -1}
  float tmp_cube_vertecies[6][4][3] = {
    {VERTEX_6, VERTEX_7, VERTEX_3, VERTEX_2}, // right
    {VERTEX_4, VERTEX_5, VERTEX_1, VERTEX_0}, // left
    {VERTEX_4, VERTEX_5, VERTEX_6, VERTEX_7}, // top
    {VERTEX_0, VERTEX_1, VERTEX_2, VERTEX_3}, // bottom
    {VERTEX_4, VERTEX_7, VERTEX_3, VERTEX_0},  // back
    {VERTEX_5, VERTEX_6, VERTEX_2, VERTEX_1} // front
  };
#undef VERTEX_0
#undef VERTEX_1
#undef VERTEX_2
#undef VERTEX_3
#undef VERTEX_4
#undef VERTEX_5
#undef VERTEX_6
#undef VERTEX_7

#define LOWER_LEFT {0, 0}
#define LOWER_RIGHT {1, 0}
#define UPPER_LEFT {0, 1}
#define UPPER_RIGHT {1, 1}
  memcpy(cube_vertecies_, tmp_cube_vertecies, sizeof(float) * 6 * 4 * 3);
  float tmp_tex_coordinate[6][4][2] = {
    {UPPER_RIGHT, UPPER_LEFT,  LOWER_LEFT,  LOWER_RIGHT},  // right
    {UPPER_RIGHT, UPPER_LEFT,  LOWER_LEFT,  LOWER_RIGHT}, // left
    {LOWER_LEFT,  UPPER_LEFT,  UPPER_RIGHT, LOWER_RIGHT}, // top
    {UPPER_LEFT,  LOWER_LEFT,  LOWER_RIGHT, UPPER_RIGHT}, // bottom
    {UPPER_LEFT,  UPPER_RIGHT, LOWER_RIGHT, LOWER_LEFT},  // back
    {UPPER_RIGHT, UPPER_LEFT,  LOWER_LEFT,  LOWER_RIGHT}  // front
  };
  memcpy(cube_tex_coordinate_, tmp_tex_coordinate, sizeof(float) * 6 * 4 * 2);
  //    for (int i = 0; i < 6; ++i) {
  //      for (int j = 0; j < 4; ++j) {
  //        tex_coordinate_[i][j][0] = (j == 0 || j == 3);
  //        tex_coordinate_[i][j][1] = (j == 0 || j == 1);
  //      }
  //    }


}

// ENVIROMENT_CUBE_MAP_H_
//    cube_vertecies_[0][0] = -0.5;
//    cube_vertecies_[0][1] = -0.5;
//    cube_vertecies_[0][2] = -0.5;
//    cube_vertecies_[1][0] = -0.5;
//    cube_vertecies_[1][1] = -0.5;
//    cube_vertecies_[1][2] = +0.5;
//    cube_vertecies_[2][0] = +0.5;
//    cube_vertecies_[2][1] = -0.5;
//    cube_vertecies_[2][2] = +0.5;
//    cube_vertecies_[3][0] = +0.5;
//    cube_vertecies_[3][1] = -0.5;
//    cube_vertecies_[3][2] = -0.5;
//    for (int i = 0; i < 4; ++i) {
//      cube_vertecies_[i + 4][0] = cube_vertecies_[i][0];
//      cube_vertecies_[i + 4][1] = cube_vertecies_[i][1] + 1.0f;
//      cube_vertecies_[i + 4][2] = cube_vertecies_[i][2];
//    }
//    float tmp_tex_coordinate[4 * 2] = {0, 0, 1.0, 0, 1, 1, 0, 1};
//    memcpy(tex_coordinate_, tmp_tex_coordinate, sizeof(float) * 4 * 2);
//    int tmp_triangle_index[12 * 3] = {
//      0, 1, 2, // bottom
//      1, 2, 3,
//      1, 6, 5, // front
//      1, 6, 3,
//      3, 7, 2, // right
//      3, 7, 6,
//      1, 4, 0, // left
//      1, 4, 5,
//      0, 7, 2, // back
//      0, 7, 4,
//      4, 6, 7, // top
//      4, 6, 5,
//    };
//    memcpy(triangle_index_, tmp_triangle_index, sizeof(int) * 12 * 3);

template class EnvironmentCubeMap<float>;
template class EnvironmentCubeMap<double>;
