#ifndef PROFILER_H
#define PROFILER_H
#include <iostream>
#include <sstream>
#include <queue>
#include <iomanip>
#include <cassert>
#include <unordered_map>
#include <stack>
#include <vector>
#include <algorithm>
#include "timer.h"

class Profiler {
public:
  typedef unsigned long long NanoSecond;
  static const NanoSecond kNullTime;
  enum {
    kMaxProfileID = 256
  };
  typedef std::unordered_map<std::string, int> Str2IntMap;
  Profiler();

  ~Profiler() {
    PrintProfileInfo(std::cout);
  }

  inline int GetId(const char* name) {
    std::string str_name(name);
    Str2IntMap::iterator pos = name2id_.find(str_name);
    if (pos == name2id_.end()) {
      name2id_[str_name] = (int) profile_names_.size();
      profile_names_.push_back(str_name);
      last_start_time_in_ns_.push_back(NanoSecond(0));
      accumulated_time_in_ns_.push_back(NanoSecond(0));
      number_of_calls_.push_back(0);
      return int(profile_names_.size()) - 1;
    } else {
      return pos->second;
    }
  }

  inline void Start(int id) {
    if (call_stack_.empty()) {
      is_root_[id] = true;
    }
    call_stack_.push(id);
    last_start_time_in_ns_[id] = world_time_->NanoSecond();
  }

  inline void Start(const char* prof_text) {
    int id = GetId(prof_text);
    Start(id);
  }

  /// prof_text not used
  inline void End(const char* prof_text = nullptr) {
    (void) prof_text;
    //    int id = GetId(prof_text);
    End(0);
  }

  /// id not used
  inline float End(int id) {
    (void) id;
    if (call_stack_.empty()) {
      std::cerr << "Profiler::EndTimer() => " << "EndTimer() called without StartTimer() called" << std::endl;
      int* err = nullptr;
      *err = 0;
      exit(0);
    }
    id = call_stack_.top();
    NanoSecond elapsed_time = (world_time_->NanoSecond() - last_start_time_in_ns_[id]);
    accumulated_time_in_ns_[id] += elapsed_time;
    number_of_calls_[id]++;
    call_stack_.pop();
    if (!call_stack_.empty()) {
      int top_id = call_stack_.top();
      if (call_graph_[top_id][id] == kNullTime) {
        call_graph_[top_id][id] = elapsed_time;
      } else {
        call_graph_[top_id][id] += elapsed_time;
      }
    }
    return elapsed_time * 1e-9f;
  }


  int GetNumOfCalls(int id);
  float GetAccumulatedTime(int id) {
    return accumulated_time_in_ns_[id] * 1e-9f;
  }
  float GetTimePerCall(int id);

  void PrintNodeInfo(float total_time, float parent_time, std::ostream& out, int id, int level);
  void PrintProfileInfo(std::ostream& out = std::cout);
private:
  WorldTime* world_time_;
  Str2IntMap name2id_;
  std::stack<int> call_stack_;
  std::vector<int> is_root_;
  NanoSecond call_graph_[kMaxProfileID][kMaxProfileID];
  std::vector<std::string> profile_names_;
  std::vector<NanoSecond> last_start_time_in_ns_;
  std::vector<NanoSecond> accumulated_time_in_ns_;
  std::vector<int> number_of_calls_;
};

inline std::ostream& operator<<(std::ostream& out, Profiler& profile) {
  profile.PrintProfileInfo(out);
  return out;
}
#endif // PROFILER_H
