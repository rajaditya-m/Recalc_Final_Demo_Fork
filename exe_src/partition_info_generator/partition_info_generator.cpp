#include <string>
#include "partitioned_mesh_generator.h"
#include "global.h"
#include "string_formatter.h"

int main(void)
{
  std::string name = dj::Format("%s/vert_partition.txt", GetPartitionFolder());
  PartitionedMeshGenerator mesh(GetMeshFile(), name.c_str());
  bool generate_vega_file = false;
//  bool generate_vega_file = true;
  mesh.WritePartitionInfo(GetPartitionFolder(), generate_vega_file);
  L("finish generating partition data");
  exit(0);
  return 0;
}

