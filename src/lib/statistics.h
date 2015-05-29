#ifndef STATISTICS_H
#define STATISTICS_H
#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>
#include <cmath>
#include <sstream>

namespace dj {
template <class T, class Float> class Statistics;

template <class T, class Float = float>
Statistics<T, Float> MinMaxSumAvgVariance(T* val, int num, int stride = 1) {
  T min = std::numeric_limits<T>::max();
  T max = std::numeric_limits<T>::min();
  T sum = T(0);
  for (int i = 0; i < num; ++i, val += stride) {
    min = std::min(min, *val);
    max = std::max(min, *val);
    sum += *val;
  }
  Float avg = sum / num;
  Float variance = 0;
  val--;
  for (int i = 0; i < num; ++i, val -= stride) {
    Float diff = *val - avg;
    variance += diff * diff;
  }
  variance /= num;
  return Statistics<T, Float>(min, max, sum, avg, variance);
}



template <class T, class Float = float>
struct Statistics {
  T min, max, sum;
  Float avg, variance;
  Statistics() {}
  Statistics(T _min, T _max, T _sum, Float _avg, Float _variance)  {
    Construct(_min, _max, _sum, _avg, _variance);
  }

  Statistics(const Statistics& other) {
    Construct(other.min, other.max, other.sum, other.avg, other.variance);
  }

  Statistics& operator=(const Statistics& other) {
    Construct(other.min, other.max, other.sum, other.avg, other.variance);
    return *this;
  }

  Statistics(T* val, int num, int stride = 1) {
    *this = MinMaxSumAvgVariance<T, Float>(val, num, stride);
  }

  void Construct(T _min, T _max, T _sum, Float _avg, Float _variance) {
    min = _min;
    max = _max;
    sum = _sum;
    avg = _avg;
    variance = _variance;
  }

};

template <class T>
struct Histogram {
  typedef std::pair<float,int> Item;
  std::vector<Item> histogram;

  Histogram(T* val, int num, int partition, int stride = 1) {
    Statistics<T, float> stat(val, num);
    float interval = (stat.max - stat.min) * 1.0f / partition;
//    std::cout << interval << std::endl;
    histogram.resize(partition + 1);
    for (int i = 0; i < partition; ++i) {
      histogram[i].first = interval * i + stat.min;
      histogram[i].second = 0;
    }
    histogram.back().first = stat.max;
    histogram.back().second = 0;

    for (int i = 0; i < num; ++i, val += stride) {
      int id = int((*val - stat.min) / interval);
      if (id >= partition) id = partition - 1;
      histogram[id].second++;
    }
  }

  Item& operator[](int idx) {
    return histogram[idx];
  }
  int Size() { return int(histogram.size()); }
};


template <class T>
void PrintHisogrm(std::vector<std::pair<T, int> >& histogram, std::ostream& out = std::cout) {
  for (int i = 0; i < int(histogram.size()) - 1; ++i) {
    std::stringstream ss;
    if (i == int(histogram.size()) - 1) {
      ss << "[" << histogram[i].first << ", " << histogram[i + 1].first << "]";
    } else {
      ss << "[" << histogram[i].first << ", " << histogram[i + 1].first << ")";
    }
    out << std::setw(30) << std::right << ss.str();
    out << "  " << histogram[i].second << "\n";
  }
}

}

template <class T, class Float = float>
std::ostream & operator<<(std::ostream& out, dj::Statistics<T, Float>& stat) {
  out << std::setw(11) << std::right << "min: " << stat.min << "\n";
  out << std::setw(11) << std::right << "max: " << stat.max << "\n";
  out << std::setw(11) << std::right << "sum: " << stat.sum << "\n";
  out << std::setw(11) << std::right << "avg: " << stat.avg << "\n";
  out << std::setw(11) << std::right << "varience: " << stat.variance << "\n";
  return out;
}

template <class T, class Float = float>
std::ostream & operator<<(std::ostream& out, dj::Histogram<T>& hist) {
  dj::PrintHisogrm(hist.histogram, out);
  return out;
}


#endif // STATISTICS_H
