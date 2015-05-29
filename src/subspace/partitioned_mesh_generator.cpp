#include <unordered_set>
#include <algorithm>
#include <set>
#include <sstream>
#include <fstream>
#include <functional>
#include "tet.h"
#include "partitioned_mesh_generator.h"
#include "tet_gen_mesh_io.h"
#include "string_formatter.h"
#include "vega_tet_mesh_io.h"
#include "print_macro.h"


PartitionedMeshGenerator::PartitionedMeshGenerator(const char *mesh_file, const char *vertex_partition_file)
  : Tet(mesh_file, 0, NULL, false) {
  std::ifstream in(vertex_partition_file);
  ASSERT(in.is_open(), P(vertex_partition_file));
  int v_num;
  in >> v_num;
  ASSERT(v_num == vertex_num_, P(v_num, vertex_num_));
  std::unordered_set<int> part_id_set;
  vertex_partition_id_.resize(vertex_num_);
  for (int i = 0; i < vertex_num_; ++i) {
    in >> vertex_partition_id_[i];
    part_id_set.insert(vertex_partition_id_[i]);
  }
  in.close();

  std::vector<int> parts(part_id_set.begin(), part_id_set.end());
  std::sort(parts.begin(), parts.end());
  part_num_ = int(parts.size());
  for (int p = 0; p < int(parts.size()); ++p) {
    ASSERT(p == parts[p]);
  }
  GeneratePartitionInfo();
  // Write vega file
  {
    std::vector<double> verts(X, X + vertex_num_ * 3);
    std::vector<int> tets(tet_, tet_ + tet_number * 4);
    VegaTetMeshIO::Instance()->Write(mesh_file, verts, tets);
  }
  ASSERT(Verify());
}

PartitionedMeshGenerator::~PartitionedMeshGenerator() {
}

bool PartitionedMeshGenerator::Verify() {
  bool success = true;
  std::vector<int> vert_map(vertex_num_, 0);
  for (int p = 0; p < part_num_; ++p) {
    for (int t = 0; t < int(domain_tets_[p].size()) / 4; ++t) {
      int* local_verts = &domain_tets_[p][t * 4];
      int global_verts[4] = {
        vert_local_id2global_id_[p][local_verts[0]],
        vert_local_id2global_id_[p][local_verts[1]],
        vert_local_id2global_id_[p][local_verts[2]],
        vert_local_id2global_id_[p][local_verts[3]],
      };
      ASSERT(vertex_partition_id_[global_verts[0]] == p);
      ASSERT(vertex_partition_id_[global_verts[1]] == p);
      ASSERT(vertex_partition_id_[global_verts[2]] == p);
      ASSERT(vertex_partition_id_[global_verts[3]] == p);
      vert_map[global_verts[0]] = 1;
      vert_map[global_verts[1]] = 1;
      vert_map[global_verts[2]] = 1;
      vert_map[global_verts[3]] = 1;
    }
  }
  int sum = 0;
  for (int v = 0; v < vertex_num_; ++v) {
    sum += vert_map[v];
    ASSERT(vert_map[v] == 1, P(v));
  }
  ASSERT(success);
  success = success && (sum == vertex_num_);
  ASSERT(success);

  for (int d = 0; d < int(domains_.size()); ++d) {
    for (int t = 0; t < int(domain_tets_[d].size()) / 4; ++t) {
      int* local_verts = &domain_tets_[d][t * 4];
      int global_verts[4] = {
        vert_local_id2global_id_[d][local_verts[0]],
        vert_local_id2global_id_[d][local_verts[1]],
        vert_local_id2global_id_[d][local_verts[2]],
        vert_local_id2global_id_[d][local_verts[3]],
      };
      int global_tet = tet_local_id2global_id_[d][t];
      ASSERT(tet_partition_info_[global_tet].first == d);
      ASSERT(tet_partition_info_[global_tet].second == t);
      int* true_global_verts = tet_ + global_tet * 4;
      for (int i = 0; i < 4; ++i) {
        ASSERT(true_global_verts[i] == global_verts[i]);
        ASSERT(dj::Abs(domain_vert_pos_[d][local_verts[i] * 3 + 0] - rest_pos_[global_verts[i] * 3 + 0]) < 1e-10);
        ASSERT(dj::Abs(domain_vert_pos_[d][local_verts[i] * 3 + 1] - rest_pos_[global_verts[i] * 3 + 1]) < 1e-10);
        ASSERT(dj::Abs(domain_vert_pos_[d][local_verts[i] * 3 + 2] - rest_pos_[global_verts[i] * 3 + 2]) < 1e-10);
      }
    }
  }
  return success;
}

void PartitionedMeshGenerator::GeneratePartitionInfo()  {
  std::unordered_set<Vec4i, Vec4iHash> domain_set;
  // tet_domain_info
  tet_domain_info_ = std::vector<std::vector<int> >(tet_number, std::vector<int>());
  // domains that contains the tet
  std::vector<Vec4i> tet_vec4i(tet_number);
  for (int t = 0; t < tet_number; ++t) {
    tet_domain_info_[t].clear();
    std::set<int> incident_domains;
    int* verts = tet_ + t * 4;
    for (int i = 0; i < 4; ++i) {
      incident_domains.insert(vertex_partition_id_[verts[i]]);
    }
    tet_domain_info_[t] = std::vector<int>(incident_domains.begin(), incident_domains.end());
    std::sort(tet_domain_info_[t].begin(), tet_domain_info_[t].end());
    Vec4i domain(tet_domain_info_[t]);
    domain_set.insert(domain);
    tet_vec4i[t] = domain;
  }

  domains_.clear();
  for (int p = 0; p < part_num_; ++p) {
    Vec4i domain(p);
    domains_.push_back(domain);
  }
  for (Vec4i d : domain_set) {
    if (d.x[1] != -1) {
      domains_.push_back(d);
    }
  }

  for (int i = 0; i < int(domains_.size()); ++i) {
    domain_index_[domains_[i]] = i;
  }
  // vertex local and global mapping
  vert_local_id2global_id_ = std::vector<std::vector<int> >(domains_.size(), std::vector<int>());
  std::vector<std::pair<int, int> > tet_domain_id(tet_number) ;
  domain_tets_.resize(domains_.size());
  domain_vert_pos_.resize(domains_.size());
  for (int t = 0; t < tet_number; ++t) {
    Vec4i& info = tet_vec4i[t];
    ASSERT(domain_index_.count(info) > 0);
    int domain_id = domain_index_[info];
    tet_domain_id[t] = std::make_pair(domain_id, t);
  }
  //  for (int t : incident_tet_[6510]) {
  //    P(t);
  //    P(tet_domain_id[t].first, tet_domain_id[t].second);
  //  }
  //  P(vertex_partition_id_[6510]);
  //  P(dj::Vec4i(&domains_[8].x[0]));
  std::sort(tet_domain_id.begin(), tet_domain_id.end());
  std::vector<int> vertex_bit_map(vertex_num_, -1);
  int prev_domain = -1;
  // vertex number and tet number in each domain
  int t_num = 0;
  int v_num = 0;

  domain_vert_pos_.clear();
  domain_tets_.clear();
  vert_local_id2global_id_.clear();
  tet_local_id2global_id_.clear();
  tet_partition_info_.clear();

  domain_vert_pos_.resize(domains_.size());
  domain_tets_.resize(domains_.size());
  vert_local_id2global_id_.resize(domains_.size());
  tet_local_id2global_id_.resize(domains_.size());
  tet_partition_info_.resize(tet_number);
  //  P(incident_tet_[6510].size());
  for (int i = 0; i < int(tet_domain_id.size()); ++i) {
    int domain_id = tet_domain_id[i].first;
    int t = tet_domain_id[i].second;
    // a new domain
    if (domain_id != prev_domain) {
      vertex_bit_map = std::vector<int>(vertex_num_, -1);
      t_num = 0;
      v_num = 0;
      prev_domain = domain_id;
      domain_vert_pos_[domain_id].clear();
      tet_local_id2global_id_[domain_id].clear();
      domain_tets_[domain_id].clear();
      vert_local_id2global_id_[domain_id].clear();
    }
    int* verts =  tet_ + t * 4;
    for (int i = 0; i < 4; ++i) {
      // new vertex in the domain
      if (vertex_bit_map[verts[i]] == -1) {
        vertex_bit_map[verts[i]] = v_num;
        vert_local_id2global_id_[domain_id].push_back(verts[i]);

        //        if (verts[i] == 6510) { P(domain_id, vert_local_id2global_id_[domain_id].size()); P(vert_local_id2global_id_[8][104]); }
        domain_vert_pos_[domain_id].push_back(rest_pos_[verts[i] * 3 + 0]);
        domain_vert_pos_[domain_id].push_back(rest_pos_[verts[i] * 3 + 1]);
        domain_vert_pos_[domain_id].push_back(rest_pos_[verts[i] * 3 + 2]);
        v_num++;
      }
      domain_tets_[domain_id].push_back(vertex_bit_map[verts[i]]);
    }
    tet_partition_info_[t] = std::make_pair(domain_id, t_num);
    tet_local_id2global_id_[domain_id].push_back(t);
    t_num++;
  }
  //  P(vert_local_id2global_id_[8][104]);

  vertex_local_id_ = std::vector<int>(vertex_num_, -1);
  for (int p = 0; p < part_num_; ++p) {
    //    P(p, vert_local_id2global_id_[p].size());
    for (int local_v = 0; local_v < int(vert_local_id2global_id_[p].size()); ++local_v) {
      int global_v = vert_local_id2global_id_[p][local_v];
      //      ASSERT(global_v != 6510, P(vert_local_id2global_id_[p][local_v], p, local_v));
      ASSERT(vertex_local_id_[global_v] == -1);
      vertex_local_id_[global_v] = local_v;
    }
  }
  for (int v = 0; v < int(vertex_local_id_.size()); ++v) {
    if (vertex_local_id_[v] == -1) {
      for (int t : incident_tet_[v]) {
        KK;
        int* verts = tet_ + t * 4;
        int vert_partition[4] = {
          vertex_partition_id_[verts[0]],
          vertex_partition_id_[verts[1]],
          vertex_partition_id_[verts[2]],
          vertex_partition_id_[verts[3]],
        };
        for (int domain : tet_domain_info_[t]) {
          P(t, domain, dj::Vec4i(tet_ + t * 4), dj::Vec4i(vert_partition));
        }
      }
      P(v, vertex_partition_id_[v]);
      ASSERT(false);
    }
  }
}

void PartitionedMeshGenerator::WritePartitionInfo(const char *folder, bool generate_vega_file) {
  auto ToString = [&](Vec4i & a) -> std::string {
    std::stringstream str;
    for (int i = 0; i < 4; ++i) {
      str << "_" << a.x[i];
      //      if (a.x[i] == -1) { str << "_-1"; } else { str << "_" << a.x[i]; }
    }
    return str.str();
  };
  auto GetFileNamePrefix = [&](int p) -> const char* {
    static char prefix[512];
    if (p < part_num_) {
      sprintf(prefix, "%s/partition_%d", folder, p);
    } else {
      sprintf(prefix, "%s/interface%s", folder, ToString(domains_[p]).c_str());
    }
    return prefix;
  };

  char file_name[512];
  // partition ids
  {
    sprintf(file_name, "%s/interface_domains.txt", folder);
    std::ofstream out(file_name);
    assert(out.is_open());
    int interface_domain_num = int(domains_.size()) - part_num_;
    out << part_num_ << "\n";
    out << interface_domain_num << "\n";
    for (int i = 0; i < interface_domain_num; ++i) {
      out << i + part_num_ << " "
          << domains_[i + part_num_].x[0] << " "
          << domains_[i + part_num_].x[1] << " "
          << domains_[i + part_num_].x[2] << " "
          << domains_[i + part_num_].x[3] << "\n";
    }
    out.close();
  }
  // tet meshes
  {
    for (int p = 0; p < int(domains_.size()); ++p) {
      sprintf(file_name, "%s", GetFileNamePrefix(p));
      TetGenMeshIO::Instance()->Write(file_name, domain_vert_pos_[p], domain_tets_[p]);
      if (generate_vega_file) {
        VegaTetMeshIO::Instance()->Write(file_name, domain_vert_pos_[p], domain_tets_[p]);
      }
    }
  }
  // vert global ids
  {
    //    for (int p = 0; p < part_num_; ++p) {
    for (int p = 0; p < int(domains_.size()); ++p) {
      sprintf(file_name, "%s.vert_global_id.txt", GetFileNamePrefix(p));
      std::ofstream out(file_name);
      assert(out.is_open());
      int v_num = int(vert_local_id2global_id_[p].size());
      out << v_num << "\n";
      for (int v = 0; v < v_num; ++v) {
        out << vert_local_id2global_id_[p][v] << "\n";
      }
      out.close();
    }
  }
  // tet global ids
  {
    //    for (int p = 0; p < part_num_; ++p) {
    for (int p = 0; p < int(domains_.size()); ++p) {
      sprintf(file_name, "%s.tet_global_id.txt", GetFileNamePrefix(p));
      std::ofstream out(file_name);
      assert(out.is_open());
      //      ASSERT(out.is_open(), P(file_name));
      int t_num = int(tet_local_id2global_id_[p].size());
      out << t_num << "\n";
      for (int t = 0; t < t_num; ++t) {
        out << tet_local_id2global_id_[p][t] << "\n";
      }
      out.close();
    }
  }
  // partition vertex mass
  {
    for (int p = 0; p < int(domains_.size()); ++p) {
      sprintf(file_name, "%s.vertex_mass.txt", GetFileNamePrefix(p));
      std::ofstream out_txt(file_name);
      sprintf(file_name, "%s.vertex_mass.bin", GetFileNamePrefix(p));
      std::ofstream out_bin(file_name, std::ios::binary);
      assert(out_txt.is_open() && out_bin.is_open());
      //      ASSERT(out.is_open(), P(file_name));
      int v_num = int(vert_local_id2global_id_[p].size());
      out_bin.write((char*) &v_num, sizeof(int));
      out_bin.write((char*) &vert_local_id2global_id_[p][0], sizeof(int) * v_num);
      out_txt << v_num << "\n";
      for (int local_v = 0; local_v < v_num; ++local_v) {
        int global_v = vert_local_id2global_id_[p][local_v];
        out_txt << mass_[global_v] << "\n";
      }
      out_txt.close();
      out_bin.close();
    }
  }
  // vert paritition info (partition_id, local vertex id)
  {
    sprintf(file_name, "%s/vertex_partition_info.txt", folder);
    std::ofstream out(file_name);
    //    ASSERT(out.is_open());
    out << vertex_num_ << "\n";
    for (int v = 0; v < vertex_num_; ++v) {
      out << vertex_partition_id_[v] << " " << vertex_local_id_[v] << "\n";
    }
    out.close();
  }
  // tet partition info (partition_id, local tet id)
  {
    sprintf(file_name, "%s/tet_partition_info.txt", folder);
    std::ofstream out(file_name);
    assert(out.is_open());
    out << tet_number << "\n";
    for (int t = 0; t < tet_number; ++t) {
      out << tet_partition_info_[t].first << " " << tet_partition_info_[t].second << "\n";
    }
    out.close();
  }
}
