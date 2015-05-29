#include <iostream>
#include <fstream>
#include <cmath>
#include "opengl_helper.h"
#include "string_formatter.h"
#include "vector_lib.h"
#include "print_macro.h"

#if 0
void DisplayString(const char* str) {
  glColor3f(0, 1, 0) ;
  glMatrixMode(GL_PROJECTION) ;
  glPushMatrix();
  glLoadIdentity();
  int screen_width = glutGet(GLUT_WINDOW_WIDTH);
  int screen_height = glutGet(GLUT_WINDOW_HEIGHT);
  gluOrtho2D(0, screen_width, 0, screen_height) ;
  glMatrixMode(GL_MODELVIEW) ;
  glPushMatrix();
  glLoadIdentity() ;
  glRasterPos2i(5, 5) ;
  for (int i = 0; str[i] != '\0'; i++) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]) ;
  }
  glMatrixMode(GL_PROJECTION) ;
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW) ;
  glPopMatrix();
}
#endif

int WritePPM(const char* output_file_name, unsigned char* image, int resolution_x, int resolution_y) {
  using namespace std;
  ofstream image_file;
  image_file.open(output_file_name) ;
  if (!image_file.is_open()) {
    std::cerr << "failed to open file for writing image." << std::endl ;
    return -1 ;
  }
  image_file << "P3" << std::endl ;
  image_file << resolution_x << " " << resolution_y << std::endl ;
  image_file << 255 << std::endl ;
  for (int y = 0 ; y < resolution_y ; y++) {
    for (int x = 0 ; x < resolution_x ; x++ ) {
      int pos = (y * resolution_x + x) * 3 ;
      image_file << ( (unsigned int) image[pos]) << " " << ( (unsigned int) image[pos + 1]) << " " << ((unsigned int) image[pos + 2]) << " " ;
    }
    image_file << std::endl ;
  }
  image_file.close() ;
  return  1 ;
}

double GetPointDepth(double* pos) {
  double winX, winY, winZ; // holder for world coordinates
  GLint view_port[4]; // viewport dimensions+pos
  GLdouble projectioin_matrix[16];
  GLdouble model_view_matrix[16];
  glGetDoublev (GL_MODELVIEW_MATRIX, model_view_matrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projectioin_matrix);
  glGetIntegerv(GL_VIEWPORT, view_port);
  gluProject(pos[0], pos[1], pos[2],
             model_view_matrix, projectioin_matrix, view_port,
             &winX, &winY, &winZ);
  return winZ;
}

float GetPointDepth(float* pos) {
  double tmp_pos[3] = {pos[0], pos[1], pos[2]};
  return (float) GetPointDepth(tmp_pos);
}

float GetPixelDepth(int pixel_position_x, int pixel_position_y) {
  GLint view_port[4];
  glGetIntegerv(GL_VIEWPORT, view_port);
  float depth;
  glReadPixels(pixel_position_x, view_port[3] - pixel_position_y, 1, 1,
               GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
  return depth;
}


void GetPixelWorldPosition(int pixel_position_x, int pixel_position_y, double depth, double* world_pos) {
  GLint view_port[4]; // viewport dimensions+pos
  GLdouble projectioin_matrix[16];
  GLdouble model_view_matrix[16];
  glGetDoublev (GL_MODELVIEW_MATRIX, model_view_matrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projectioin_matrix);
  glGetIntegerv(GL_VIEWPORT, view_port);
  //view[3]-cursorY = conversion from upper left (0,0) to lower left (0,0)
  //Unproject 2D Screen coordinates into wonderful world coordinates
  gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, depth,
               model_view_matrix, projectioin_matrix, view_port,
               world_pos + 0, world_pos + 1, world_pos + 2);
}

void GetPixelWorldPosition(int pixel_position_x, int pixel_position_y, float depth, float* world_pos) {
  double tmp_pos[3];
  GetPixelWorldPosition(pixel_position_x, pixel_position_y, (double) depth, tmp_pos);
  world_pos[0] = float(tmp_pos[0]);
  world_pos[1] = float(tmp_pos[1]);
  world_pos[2] = float(tmp_pos[2]);
}

void GetSelectionRay(int pixel_position_x, int pixel_position_y, double *starting_point, double *ending_point) {
  double objX, objY, objZ; // holder for world coordinates
  GLint view_port[4]; // viewport dimensions+pos
  GLdouble projectioin_matrix[16];
  GLdouble model_view_matrix[16];
  glGetDoublev (GL_MODELVIEW_MATRIX, model_view_matrix);
  glGetDoublev (GL_PROJECTION_MATRIX, projectioin_matrix);
  glGetIntegerv(GL_VIEWPORT, view_port);
  //view[3]-cursorY = conversion from upper left (0,0) to lower left (0,0)
  //Unproject 2D Screen coordinates into wonderful world coordinates
  gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, 1,
               model_view_matrix, projectioin_matrix, view_port, &objX, &objY, &objZ);
  ending_point[0] = objX;
  ending_point[1] = objY;
  ending_point[2] = objZ;
  gluUnProject((GLdouble) pixel_position_x, (GLdouble) view_port[3] - pixel_position_y, 0,
               model_view_matrix, projectioin_matrix, view_port, &objX, &objY, &objZ);
  starting_point[0] = objX;
  starting_point[1] = objY;
  starting_point[2] = objZ;
}

void GetSelectionRay(int pixel_position_x, int pixel_position_y, float *starting_point, float *ending_point) {
  double start[3], end[3];
  GetSelectionRay(pixel_position_x, pixel_position_y, start, end);
  starting_point[0] = (float) start[0];
  starting_point[1] = (float) start[1];
  starting_point[2] = (float) start[2];
  ending_point[0] = (float) end[0];
  ending_point[1] = (float) end[1];
  ending_point[2] = (float) end[2];
}

void IntersectWithHorizontalPlane(int pixel_position_x, int pixel_position_y, float &pos_x, float &pos_y) {
  float starting_point[3];
  float ending_point[3];
  GetSelectionRay(pixel_position_x, pixel_position_y, starting_point, ending_point);
  pos_x = -starting_point[1] / (ending_point[1] - starting_point[1]) * (ending_point[0] - starting_point[0]) + starting_point[0];
  pos_y = -starting_point[1] / (ending_point[1] - starting_point[1]) * (ending_point[2] - starting_point[2]) + starting_point[2];
}

void DrawSphere(float radius, int lats, int longs) {

  GLUquadricObj* Sphere = gluNewQuadric();
  gluSphere(Sphere, radius, lats, longs);
  //  gluDeleteQuadric(Sphere);
  //  const float kPi = 3.1415926f;
  //  for (int i = 0; i <= lats; i++) {
  //    float lat0 = kPi * (-0.5 + (float) (i - 1) / lats);
  //    float z0  = radius * sin(lat0);
  //    float zr0 =  radius* cos(lat0);
  //
  //    float lat1 = kPi * (-0.5 + (float) i / lats);
  //    float z1 = radius * sin(lat1);
  //    float zr1 = radius * cos(lat1);
  //
  //    glBegin(GL_QUAD_STRIP);
  //    for (int j = 0; j <= longs; j++) {
  //      float lng = 2 * kPi * (float) (j - 1) / longs;
  //      float x = cos(lng);
  //      float y = sin(lng);
  //      glNormal3f(x * zr0, y * zr0, z0);
  //      glVertex3f(x * zr0, y * zr0, z0);
  //      glNormal3f(x * zr1, y * zr1, z1);
  //      glVertex3f(x * zr1, y * zr1, z1);
  //    }
  //    glEnd();
  //  }
}


void DrawCylinder(float radius, float height, int slice, int stack) {
  static GLUquadricObj *quadObj = gluNewQuadric();
  gluCylinder(quadObj,
              radius, // base radius
              radius, // top radius
              height, // height
              slice, // slice
              stack // stack
             );
}


void DrawGradientBackGround(float *lower_color, float *upper_color) {
  // Maya background
  float default_lower_color[3] = {30 / 255.0f, 30 / 255.0f, 30 / 255.0f};
  float default_upper_color[3] = {130 / 255.0f, 140 / 255.0f, 160 / 255.0f};
  // Meshlab backgroud
  //  float default_lower_color[3] = {115 / 255.0f, 115 / 255.0f, 230 / 255.0f};
  //  float default_upper_color[3] = {0 / 255.0f, 0 / 255.0f, 0 / 255.0f};
  // White back ground
  //  float default_lower_color[3] = {250 / 255.0f, 250 / 255.0f, 250 / 255.0f};
  //  float default_upper_color[3] = {220 / 255.0f, 220 / 255.0f, 220 / 255.0f};

  if (lower_color == NULL) { lower_color = default_lower_color; }
  if (upper_color == NULL) { upper_color = default_upper_color; }
  glPushAttrib(GL_ENABLE_BIT);
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);
  glBegin(GL_QUADS);
  glColor3fv(lower_color);
  glVertex3f(-1.0, -1.0, 0);
  glVertex3f(1.0, -1.0, 0);
  glColor3fv(upper_color);
  glVertex3f(1.0, 1.0, 0);
  glVertex3f(-1.0, 1.0, 0);
  glEnd();
  glPopAttrib();

  glPopMatrix(); // Pop model view matrix
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
}


void DrawAxis() {
  GLboolean lighting_enabled;
  glGetBooleanv(GL_LIGHTING, &lighting_enabled);
  glDisable(GL_LIGHTING);
  // X axis
  glColor3f(1, 0, 0) ;
  DrawArrow<float>(0.0, 0, 0, 1, 0, 0);
  //  // Y axis
  glColor3f(0, 1, 0) ;
  DrawArrow<float>(0.0, 0, 0, 0, 1, 0);
  //  // Z axis
  glColor3f(0, 0, 1);
  DrawArrow<float>(0.0, 0, 0, 0, 0, 1);
  if (lighting_enabled) { glEnable(GL_LIGHTING);}
}

///////////////////////////////////////////////////////////////////////////////
// write 2d text using GLUT
// The projection matrix must be set to orthogonal before call this function.
///////////////////////////////////////////////////////////////////////////////
//void drawString(const char *str, int x, int y, float color[4], void *font)
//{
//    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT); // lighting and color mask
//    glDisable(GL_LIGHTING);     // need to disable lighting for proper text color
//    glDisable(GL_TEXTURE_2D);

//    glColor4fv(color);          // set text color
//    glRasterPos2i(x, y);        // place text position

//    // loop all characters in the string
//    while(*str)
//    {
//        glutBitmapCharacter(font, *str);
//        ++str;
//    }

//    glEnable(GL_TEXTURE_2D);
//    glEnable(GL_LIGHTING);
//    glPopAttrib();
//}




void DrawCheckBoard(float width, float height,
                    int x_slice, int y_slice,
                    const float *color1, const float *color2,
                    bool enable_lighting) {
  float x_slice_size = width / x_slice;
  float y_slice_size = height / y_slice;
  glPushAttrib(GL_ENABLE_BIT);
  if (enable_lighting) {
    glEnable(GL_LIGHTING);
  }
  glEnable(GL_DEPTH_TEST);
  glBegin(GL_QUADS);
  float x = -width / 2;// + i * x_slice_size;
  for (int i = 0; i < x_slice + 1; ++i, x += x_slice_size) {
    float y = -height / 2;// + j * y_slice_size;
    for (int j = 0; j < y_slice + 1; ++j, y += y_slice_size) {
      const float* color = ((i + j) % 2 == 0) ? color1 : color2;
      if (enable_lighting) {
        float diffuse[4] = {color[0], color[1], color[2], 1};
        //        float diffuse[4] = {0.4, 0.8, 0.2, 1};
        float ambient[4] = {0.2f, 0.2f, 0.2f, 1};
        //        float specular[4] = {1, 1, 1, 1};
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
        //        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
        glNormal3f(0, 0, 1);
      } else {
        glColor3fv(color);
      }
      glVertex2f(x, y);
      glVertex2f(x + x_slice_size, y);
      glVertex2f(x + x_slice_size, y + y_slice_size);
      glVertex2f(x, y + y_slice_size);
    }
  }
  glEnd();
  glPopAttrib();
}


//void DrawString3D(const char *str, const float pos[], const float color[], void *font)
//{
//  glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT); // lighting and color mask
//  glDisable(GL_LIGHTING);     // need to disable lighting for proper text color
//  glDisable(GL_TEXTURE_2D);

//  glColor3fv(color);          // set text color
//  glRasterPos3fv(pos);        // place text position

//  // loop all characters in the string
//  while(*str)
//  {
//    glutBitmapCharacter(font, *str);
//    ++str;
//  }

//  glEnable(GL_TEXTURE_2D);
//  glEnable(GL_LIGHTING);
//  glPopAttrib();
//}


void GetPixelWorldPosition(int pixel_position_x, int pixel_position_y, double *world_pos) {
  double depth = GetPixelDepth(pixel_position_x, pixel_position_y);
  GetPixelWorldPosition(pixel_position_x, pixel_position_y, depth, world_pos);
}


void GetPixelWorldPosition(int pixel_position_x, int pixel_position_y, float *world_pos) {
  float depth = GetPixelDepth(pixel_position_x, pixel_position_y);
  GetPixelWorldPosition(pixel_position_x, pixel_position_y, depth, world_pos);
}


void Texture2NormalMap(unsigned char *texture, int width, int height, float scale) {
  auto GetIntensity = [](unsigned char * pixel) {
    const float r = float(pixel[0]);
    const float g = float(pixel[1]);
    const float b = float(pixel[2]);
    const float average = (r + g + b) / 3.0f;
    return average / 255.0f;
  };
  auto Clamp = [](int x, int max) {
    if (x > max) return max;
    else if (x < 0) return 0;
    else return x;
  };
  // transform -1 - 1 to 0 - 255
  auto ToRgb = [](float * x, unsigned char * pixel) {
    pixel[0] = (unsigned char)((x[0] + 1.0f) * (255.0f / 2.0f));
    pixel[1] = (unsigned char)((x[1] + 1.0f) * (255.0f / 2.0f));
    pixel[2] = (unsigned char)((x[2] + 1.0f) * (255.0f / 2.0f));
  };
  auto Idx = [&](int x, int y) {
    x = Clamp(x, width - 1);
    y = Clamp(y, height - 1);
    return ((y * width) + x);
  };

  std::vector<float> intensity_array(width * height);
  for (int y = 0, idx = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x, ++idx) {
      intensity_array[idx] = GetIntensity(texture + idx * 3);
    }
  }

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      float intensity[3][3];
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          intensity[i][j] = intensity_array[Idx(x + (i - 1), y + (j - 1))];
        }
      }
      // Sobel filter
      float dxyz[3] = {
        (intensity[2][2] + 2.0f * intensity[1][2] + intensity[0][2]) - (intensity[2][0] + 2.0f * intensity[1][0] + intensity[0][0]),
        (intensity[0][0] + 2.0f * intensity[0][1] + intensity[0][2]) - (intensity[2][0] + 2.0f * intensity[2][1] + intensity[2][2]),
        1.0f / scale,
      };
      dj::Normalize3(dxyz);
      ToRgb(dxyz, texture + Idx(x, y) * 3);
    }
  }
}

//**************************************************************************************
//  Read shaders from the disk into the main memory
//**************************************************************************************
int Read_Shader(const char *name, char **shader_text) {
  FILE* fp = fopen(name, "r+");
  if (fp == NULL) return 0;
  //Calculate the file size and allocate the shader_content
  fseek(fp, 0L, SEEK_END);
  int size = ftell(fp) + 1;
  if (*shader_text)	delete[] (*shader_text);
  *shader_text = new char[size];
  //Read the shader file
  int count = 0;
  fseek(fp, 0, SEEK_SET);
  count = (int) fread(*shader_text, 1, size, fp);
  (*shader_text)[count] = '\0';
  //The end
  fclose(fp);
  return count;
}

bool Read_Shader_Source(const char *shader_name, GLchar **vertexShader, GLchar **fragmentShader) {
  char vert_shader_name[1024];
  sprintf(vert_shader_name, "%s.vert", shader_name);
  if (!Read_Shader(vert_shader_name, vertexShader)) {
    printf("Cannot read the file %s\n", vert_shader_name);
    return false;
  }
  char frag_shader_name[1024];
  sprintf(frag_shader_name, "%s.frag", shader_name);
  if (!Read_Shader(frag_shader_name, fragmentShader)) {
    printf("Cannot read the file %s\n", frag_shader_name);
    return false;
  }
  return true;
}


//**************************************************************************************
//  GLSL setup routine
//**************************************************************************************
GLuint SetupGLSL(const char *shader_file_name) {
  //Step 2: Create the objects
  GLuint programObject;
  GLuint vertexShaderObject;
  GLuint fragmentShaderObject;
  if (!(programObject = glCreateProgram())) {
    printf("Error creating shader program object.\n");
    exit(1);
  } else {
    //    printf("Succeeded creating shader program object.\n");
  }
  if (!(vertexShaderObject = glCreateShader(GL_VERTEX_SHADER))) {
    printf("Error creating vertex shader object.\n");
    exit(1);
  } else {
    //    printf("Succeeded creating vertex shader object.\n");
  }
  if (!(fragmentShaderObject = glCreateShader(GL_FRAGMENT_SHADER))) {
    printf("Error creating fragment shader object.\n");
    exit(1);
  } else {
    //    printf("Succeeded creating fragment shader object.\n");
  }

  //Step 3: Load the shaders from the disk into the two shader objects
  GLchar* vertexShaderSource = 0;
  GLchar* fragmentShaderSource = 0;
  Read_Shader_Source(shader_file_name, &vertexShaderSource, &fragmentShaderSource);
  glShaderSource(vertexShaderObject, 1, (const GLchar**)&vertexShaderSource, NULL);
  glShaderSource(fragmentShaderObject, 1, (const GLchar**)&fragmentShaderSource, NULL);
  delete[] vertexShaderSource;
  delete[] fragmentShaderSource;

  //Step 4: Compile the vertex shader
  glCompileShader(vertexShaderObject);
  //If there is any error, print out the error log
  GLint result;
  glGetShaderiv(vertexShaderObject, GL_COMPILE_STATUS, &result);
  if (result == GL_FALSE) {
    printf(" vertex shader compilation failed!\n");
    GLint logLen;
    glGetShaderiv(vertexShaderObject, GL_INFO_LOG_LENGTH, &logLen);
    if (logLen > 0) {
      char* log = new char[logLen];
      GLsizei written;
      glGetShaderInfoLog(vertexShaderObject, logLen, &written, log);
      printf("Shader log: \n %s", log);
      delete[] log;
    }
    auto msg = dj::Format("Failed to compile vertex shaeder \"%z\"", shader_file_name);
    ASSERT(false, L(msg));
  }

  //Step 5: Compile the fragment shader
  glCompileShader(fragmentShaderObject);
  //If there is any error, print out the error log
  glGetShaderiv(fragmentShaderObject, GL_COMPILE_STATUS, &result);
  if (result == GL_FALSE) {
    printf(" fragment shader compilation failed!\n");
    GLint logLen;
    glGetShaderiv(fragmentShaderObject, GL_INFO_LOG_LENGTH, &logLen);
    if (logLen > 0) {
      char* log = new char[logLen];
      GLsizei written;
      glGetShaderInfoLog(fragmentShaderObject, logLen, &written, log);
      printf("Shader log: \n %s", log);
      delete[] log;
    }
    auto msg = dj::Format("Failed to compile fragment shaeder \"%z\"", shader_file_name);
    ASSERT(false, L(msg));
  }

  //Step 6: Attach the shader objects and link the program object
  glAttachShader(programObject, vertexShaderObject);
  glAttachShader(programObject, fragmentShaderObject);
  glLinkProgram(programObject);

  return (programObject);
}

int CreateTexture(unsigned char* image, int width, int height,
                  GLenum mag_filter, GLenum min_filter, GLint pixel_format,
                  GLenum wrap_s, GLenum wrap_t) {
  GLuint texture_id;
  glGenTextures(1, &texture_id);
  glBindTexture(GL_TEXTURE_2D, texture_id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap_s);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap_t);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
  glTexImage2D(GL_TEXTURE_2D, 0, pixel_format, width, height, 0, pixel_format, GL_UNSIGNED_BYTE, image);
  glBindTexture(GL_TEXTURE_2D, 0);
  return texture_id;
}


/// error code strings come from https://www.opengl.org/wiki/OpenGL_Error
const char *GetGLErrorString() {
  int code = glGetError();
  switch (code) {
    case GL_NO_ERROR:
      return "NO ERROR";
    case 0x0500:
      return
        "GL_INVALID_ENUM​, 0x500,\n"
        "Given when an enumeration parameter is not a legal enumeration for that function."
        "This is given only for local problems; if the spec allows the enumeration in certain circumstances,"
        "where other parameters or state dictate those circumstances, then GL_INVALID_OPERATION​ is the result instead.";
    case 0x0501:
      return
        "GL_INVALID_VALUE​, 0x0501\n"
        "Given when a value parameter is not a legal value for that function."
        "This is only given for local problems; if the spec allows the value in certain circumstances,"
        "where other parameters or state dictate those circumstances, then GL_INVALID_OPERATION is the result instead." ;
    case 0x0502:
      return
        "GL_INVALID_OPERATION​, 0x0502\n"
        "Given when the set of state for a command is not legal for the parameters given to that command."
        "It is also given for commands where combinations of parameters define what the legal parameters are.";
    case 0x0503:
      return
        "GL_STACK_OVERFLOW​, 0x0503\n"
        "Given when a stack pushing operation cannot be done because it would overflow the limit of that stack's size";
    case 0x0504:
      return
        "GL_STACK_UNDERFLOW​, 0x0504\n"
        "Given when a stack popping operation cannot be done because the stack is already at its lowest point.";
    case 0x0505:
      return
        "GL_OUT_OF_MEMORY​, 0x0505\n"
        "Given when performing an operation that can allocate memory, and the memory cannot be allocated."
        "The results of OpenGL functions that return this error are undefined; it is allowable for partial operations to happen.";
    case 0x0506:
      return
        "GL_INVALID_FRAMEBUFFER_OPERATION​, 0x0506\n"
        "Given when doing anything that would attempt to read from or write/render to a framebuffer that is not complete.";
    case 0x0507:
      return
        "GL_CONTEXT_LOST​, 0x0507 (with OpenGL 4.5 or ARB_KHR_robustness)\n"
        "Given if the OpenGL context has been lost, due to a graphics card reset.";
    case 0x8031:
      return
        "GL_TABLE_TOO_LARGE​1, 0x8031\n"
        "Part of the ARB_imaging extension.";
    default:
      return "Invalid error coe";
  }
}




