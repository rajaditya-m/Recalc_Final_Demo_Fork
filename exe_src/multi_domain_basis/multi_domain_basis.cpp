#include <string>
#include "print_macro.h"

#ifdef HAS_BASIS_GENERATION_MODULE
#include "string_formatter.h"
#include "config_file.h"
#include "global.h"
#include "skeleton_node.h"
#include "skeleton.h"
#include "multi_domain_basis_generator.h"
#endif

int main(void) {
#ifdef HAS_BASIS_GENERATION_MODULE
  Skeleton skeleton(global::skeleton_file.c_str(), true);
  (void) skeleton;
  std::string mesh_file = dj::Format("%z/%z", GetModelFolder(), pose_conf.Get<std::string>("model name"));
  P(mesh_file);
  //    std::string vertex_partition_file = dj::Format("%z/vertex_partition_info.txt", GetPartitionFolder());
  MultiDomainBasisGenerator gen(mesh_file.c_str(), GetPartitionFolder());
  gen.SetFixedVertex(3);
//      gen.AttachSkeleton(global::skeleton);
  auto folder = dj::Format("%z/modal_basis", GetDataFolder());
  gen.GenerateBasis(folder.c_str(), 26, 20);
//  gen.GenerateBasis(folder.c_str(), 1, 1);
  L("basis generated.");
#else
  L("Basis generation modules is not enabled");
#endif
  return 0;
}

