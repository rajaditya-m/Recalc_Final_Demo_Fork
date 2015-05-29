#include <iostream>
#include <string>
#include "multi_domain_cubature.h"
#include "single_domain_cubature.h"
#include "string_formatter.h"
#include "global.h"
#include "print_macro.h"

int main(void) {
    float scale = 0.01000;
    int maxCubaturePoints = 60;
    SingleDomainCubature singleDomainCubOp(GetMeshFile());
    std::string output_folder = dj::Format("%z/cubature", GetDataFolder());
    std::string name = dj::Format("%s/modal_basis/basis_100.U", GetDataFolder());
    singleDomainCubOp.LoadBinarySubspace(name.c_str());
    singleDomainCubOp.SetFolder(100,output_folder,scale,true);
    singleDomainCubOp.GenerateCubature(maxCubaturePoints, 1.0e-6);
}

