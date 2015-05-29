#include "global.h"
#include "macro_constant.h"
#include "config_file.h"
#include <string>
#include <cstdio>
#include "string_formatter.h"
class Tet;
class Skeleton;
Profiler profiler;
ConfigFile conf(DATA_DIRECTORY "../conf/armadillo.conf");
//ConfigFile conf(DATA_DIRECTORY "armadillo.conf");

namespace global
{
#if defined(_WIN32) || defined(_WIN64)
const char kTmpFolder[] = "C:/tmp/render";
#elif defined(__APPLE__)
const char kTmpFolder[] = "/Volumes/ram";
#else
const char kTmpFolder[] = "/tmp/log";
#endif
OpenGLQt* gl = NULL;
const char kPWD[] = DATA_DIRECTORY;
const char kDataDir[] = DATA_DIRECTORY;
std::string data_directory(DATA_DIRECTORY);

int simulation_step_per_idle = conf.Get<int>("simulation_step_per_idle");
int simulate = conf.Get<int>("simulate");
int pause_per_frame = conf.Get<int>("pause_per_frame");
int sim_state = 1;
float ap1 = 0.0018;
float ap2 = 0.0005;
int stiffyMult = 1;
Real time_step = conf.Get<Real>("time_step");
Real* gravity = conf.Get<Real*>("gravity");


//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// posing
Tet* current_body = NULL;

Skeleton* skeleton;
std::string skeleton_file = global::data_directory + conf.Get<std::string>("skeleton_file");
}

namespace rod
{
int pbd_iteration = conf.Get<int>("pbd_iteration");
Real* bending_twising_stiffness = conf.Get<Real*>("bending_twisting_stiffness");
int rod_v_num = conf.Get<int>("rod_v_num");
Real rod_length = conf.Get<Real>("rod_length");
}

namespace body {
const float kMinJointLength = 0.2f / 6;
const float* kBodyColor = conf.Get<float*>("kBodyColor");//{205 / 255.0f, 156 / 255.0f, 98 / 255.0f, 1};
}

const char* GetModelFolder() {
  static std::string folder = dj::GetStr("%s%s", DATA_DIRECTORY, conf.Get<std::string>("model name").c_str());
  return folder.c_str();
}

const char* GetPartitionFolder() {
  static std::string folder = dj::GetStr("%s%s/partition", DATA_DIRECTORY, conf.Get<std::string>("model name").c_str());
  return folder.c_str();
}

const char* GetDataFolder() {
  static std::string folder = dj::GetStr("%s%s/data", DATA_DIRECTORY, conf.Get<std::string>("model name").c_str());
  return folder.c_str();
}

const char* GetMeshFile() {
  static std::string file = dj::Format("%s%z/%z", DATA_DIRECTORY,
                                       conf.Get<std::string>("model name"),
                                       conf.Get<std::string>("model name"));
  return file.c_str();
}

