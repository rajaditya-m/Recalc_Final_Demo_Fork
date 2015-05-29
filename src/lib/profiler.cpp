#include "profiler.h"
#include "print_macro.h"
const Profiler::NanoSecond Profiler::kNullTime = NanoSecond(-1);

Profiler::Profiler() {
  world_time_ = WorldTime::Instance();
  accumulated_time_in_ns_.reserve(kMaxProfileID);
  number_of_calls_.reserve(kMaxProfileID);
  std::fill_n(&call_graph_[0][0], kMaxProfileID * kMaxProfileID, Profiler::kNullTime);
  is_root_ = std::vector<int>(kMaxProfileID, false);
}

int Profiler::GetNumOfCalls(int id) {
  return number_of_calls_[id];
}

float Profiler::GetTimePerCall(int id) {
  if (number_of_calls_[id] != 0) {
    return GetAccumulatedTime(id) / number_of_calls_[id];
  } else {
    std::cerr << "Profiler::GetTimePerCall() => " << profile_names_[id] << " has never been profiled" << std::endl;
    int* err = nullptr;
    *err = 0;
    exit(0);
  }
}

void Profiler::PrintNodeInfo(float total_time, float parent_time, std::ostream &out, int id, int level) {
  //    const char kLevelSymbol[] = {'-', '+', '*', '=', '#', '!', '$', '`', '~', '^', '&'};
  //    const int kNumSymbol = sizeof(kLevelSymbol) / sizeof(char);
  double percentage_total = GetAccumulatedTime(id) * 100.0 / total_time;
  double percentage_wrt_parent = GetAccumulatedTime(id) * 100.0 / parent_time;
  std::stringstream stream_total;
  stream_total << percentage_total << "%";
  std::stringstream stream_parent;
  stream_parent << percentage_wrt_parent << "%";
  out << std::setw(20) << std::left << number_of_calls_[id]
      << std::setw(20) << std::left << GetAccumulatedTime(id) / number_of_calls_[id]
      << std::setw(20) << std::left << GetAccumulatedTime(id)
      << std::setw(20) << std::left << stream_parent.str()
      << std::setw(20) << std::left << stream_total.str();
  for (int k = 0; k < level; ++k) {
    //      const char symbol = kLevelSymbol[(level - 1) % kNumSymbol];
    //      out << symbol << symbol;
    out << "  ";
  }
  out << std::left << profile_names_[id] << std::endl;
}

void Profiler::PrintProfileInfo(std::ostream &out) {
  std::vector<std::pair<double, int> > all_info;
  all_info.reserve(accumulated_time_in_ns_.size());
  NanoSecond total_time_in_ns = 0;
  for (unsigned i = 0; i < accumulated_time_in_ns_.size(); ++i) {
    if (is_root_[i]) {
      total_time_in_ns += accumulated_time_in_ns_[i];
      all_info.push_back(std::make_pair(-GetAccumulatedTime(i), int(i)));
    }
  }
  double total_time = total_time_in_ns * 1e-9f;
  std::sort(all_info.begin(), all_info.end());
  out << std::setw(20) << std::left << "num. of calls"
      << std::setw(20) << std::left << "time per call"
      << std::setw(20) << std::left << "accumulated time"
      << std::setw(20) << std::left << "% wrt. parent"
      << std::setw(20) << std::left << "% wrt. total"
      << std::setw(20) << std::left << "name"
      << std::endl;
  for (unsigned i = 0; i < all_info.size(); ++i) {
    int idx = all_info[i].second;
    std::vector<int> visited(kMaxProfileID, false);
    visited[idx] = true;
    std::stack<int> stack;
    stack.push(idx);
    PrintNodeInfo(total_time, total_time, out, idx, 0);
    int level = 0;
    while (!stack.empty()) {
      int top = stack.top();
      int max_child = -1;
      NanoSecond max_time = NanoSecond(0);
      for (int k = 0; k < kMaxProfileID; ++k) {
        if (!visited[k] && call_graph_[top][k] != kNullTime && accumulated_time_in_ns_[k] > max_time) {
          max_child = k;
          max_time = accumulated_time_in_ns_[k];
        }
      }
      if (max_child != -1) {
        stack.push(max_child);
        visited[max_child] = true;
        level++;
        PrintNodeInfo(total_time, GetAccumulatedTime(top), out, max_child, level);
      } else {
        stack.pop();
        level--;
      }
    }
  }
}
