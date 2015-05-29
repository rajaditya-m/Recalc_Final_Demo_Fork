#ifndef BASIS_IO_H_
#define BASIS_IO_H_
//#pragma once
#include <vector>

int WriteBasisInBinary(const char *file_name, int vertex_num, int basis_num, double *basis);
int WriteBasisInText(const char *file_name, int vertex_num, int basis_num, double *basis);

int ReadBasisInBinary(const char *file_name, int &vertex_num, int &basis_num, std::vector<double> &basis);
int ReadBasisInText(const char *file_name, int &vertex_num, int &basis_num, std::vector<double> &basis);

int loadCommaList(const char * filename, int * numListEntries, int ** listEntries, int offset=0);
#endif
