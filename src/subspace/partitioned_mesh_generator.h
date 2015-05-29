#ifndef PARTITIONEDMESHGENERATOR_H
#define PARTITIONEDMESHGENERATOR_H
#include <sstream>
#include <iostream>
#include <utility>
#include <unordered_map>
#include <vector>
#include "tet.h"

class PartitionedMeshGenerator : public Tet {
public:
  struct Vec4i {
    int x[4];
    Vec4i(int _x = -1, int _y = -1, int _z = -1, int _w = -1) {
      x[0] = _x; x[1] = _y; x[2] = _z; x[3] = _w;
    }
    Vec4i(const std::vector<int>& vec) {
      x[0] = x[1] = x[2] = x[3] = -1;
      for (int i = 0; i < 4 && i < int(vec.size()); ++i) {
        x[i] = vec[i];
      }
    }

    std::string ToString() const {
      std::stringstream out;
      for (int i = 0; i < 4; ++i) {
        if (x[i] == -1) {
          out << "X ";
        } else {
          out << x[i] << " ";
        }
      }
      return out.str();
    }

    bool operator==(const Vec4i& other) const {
      return  x[0] == other.x[0] && x[1] == other.x[1] && x[2] == other.x[2] && x[3] == other.x[3];
    }

    std::size_t operator()() const {
      std::hash<int> hash;
      return hash(x[0] << 24 | x[1] << 16 | x[2] << 8 | x[3]);
    }
  };

  struct Vec4iHash {
    std::size_t operator()(const Vec4i& k) const {
      return k();
    }
  };

  PartitionedMeshGenerator(const char *mesh_file, const char* vertex_partition_file);
  virtual ~PartitionedMeshGenerator();
  bool Verify();
  void GeneratePartitionInfo();
  void WritePartitionInfo(const char* folder, bool generate_vega_file = false);

  int domain_num_; // number of domains including interface domains
  int part_num_; // number of inner domains

  std::unordered_map<Vec4i, int, Vec4iHash> domain_index_;
  std::vector<Vec4i> domains_;
  //  std::unordered_map<int, int> part_id2compact_id_;
  //  std::vector<int> parts_;
  std::vector<std::vector<double> > domain_vert_pos_; // vertex posiiton for each domain (including interface domain)
  std::vector<std::vector<int> > domain_tets_; // list tets for each domain
  std::vector<std::vector<int> > vert_local_id2global_id_;
  std::vector<std::vector<int> > tet_domain_info_; // list of domain that contains the tet vertex
  std::vector<std::pair<int, int> > tet_partition_info_; // (domain_id, tet_id) for each tet
  std::vector<std::vector<int> > tet_local_id2global_id_;
  std::vector<int> vertex_local_id_; // vertex local id in each inner domain
  std::vector<int> vertex_partition_id_;
};

#endif // PARTITIONEDMESHGENERATOR_H
