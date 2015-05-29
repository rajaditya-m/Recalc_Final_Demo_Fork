#ifndef GLOBAL_H_
#define GLOBAL_H_
#include <vector>
#include <string>
#define OPEN_MP
#ifdef OPEN_MP
#  include <omp.h>
#  ifndef _MSC_VER
#    define OMP_FOR _Pragma("omp parallel for")
#  else
#    define OMP_FOR __pragma(omp parallel for)
#  endif // vs compiler
#else // no omp
#  define OMP_FOR
#endif // #ifndef OPEN_MP

#ifdef _MSC_VER
#undef _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
//#define _CRT_SECURE_NO_DEPRECATE
#endif // _MSC_VER

#if defined(_WIN32) || defined(_WIN64)
#define TOKEN_TO_STRING1(t) # t
#define STRINGIZE_TOKEN1(t) TOKEN_TO_STRING1(t)
#define DATA_DIRECTORY STRINGIZE_TOKEN1(CURRENT_DIRECTORY)
#endif
class OpenGLQt;
typedef double Real;

#include "profiler.h"
class Profiler;
extern Profiler profiler;


class Tet;
class TetSampler;
class Skeleton;
class ConfigFile;
extern ConfigFile conf;
extern ConfigFile conf;

namespace global
{
extern OpenGLQt* gl;
extern const char kTmpFolder[];
extern int simulate;
extern int sim_state;
extern int pause_per_frame;
extern const char kPWD[];
extern const char kDataDir[];
extern float ap1;
extern float ap2;
extern int stiffyMult;
extern std::string data_directory;
extern Real time_step;
extern Real* gravity;
extern int simulation_step_per_idle;


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// posing
static const float kPi = 3.1415926f;
extern Tet* current_body;
extern std::string skeleton_file;
extern Skeleton* skeleton;
}

namespace rod {
extern Real* bending_twising_stiffness;
extern int rod_v_num;
extern Real rod_length;
extern int pbd_iteration;
}


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
namespace body
{
extern const float* kBodyColor;
extern const float kMinJointLength;
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


const char* GetModelFolder();
const char* GetPartitionFolder();
const char* GetDataFolder();
const char* GetMeshFile();
namespace arch = global;
#endif // GLOBAL_H_
