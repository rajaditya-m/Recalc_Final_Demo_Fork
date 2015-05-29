#include "basis_io.h"
#include "matrixIO.h"
#include <fstream>
#include "print_macro.h"

int WriteBasisInBinary(const char *file_name, int vertex_num, int basis_num, double *basis) {
  return WriteMatrixToDisk(file_name, 3 * vertex_num, basis_num, basis);
}

int WriteBasisInText(const char *file_name, int vertex_num, int basis_num, double *basis) {
  std::ofstream out(file_name);
  ASSERT(out.is_open(), P(file_name));
  out << vertex_num * 3 << " ";
  out << basis_num << "\n";
  for (int r = 0; r < basis_num; ++r) {
    for (int v = 0; v < vertex_num * 3; ++v) {
      out << basis[r * vertex_num * 3 + v] << " ";
    }
    out << "\n";
  }
  out.close();
  return 0;
}

int ReadBasisInBinary(const char *file_name, int &vertex_num, int &basis_num, std::vector<double> &basis) {
  std::ifstream in(file_name, std::ios::binary);
  ASSERT(in.is_open(), P(file_name));
  int row_num;
  in.read((char*) &row_num, sizeof(int));
  vertex_num = row_num / 3;
  in.read((char*) &basis_num, sizeof(int));
  basis.resize(row_num * basis_num);
  in.read((char*) &basis[0], sizeof(double) * row_num * basis_num);
  in.close();
  return 0;
}

int ReadBasisInText(const char *file_name, int &vertex_num, int &basis_num, std::vector<double> &basis) {
  std::ifstream in(file_name);
  ASSERT(in.is_open(), P(file_name));
  int row_num;
  in >> row_num >> basis_num;
  vertex_num = row_num / 3;
  basis.resize(row_num * basis_num);
  for (int i = 0; i < row_num * basis_num; ++i) {
    in >> basis[i];
  }
  in.close();
  return 0;
}

// removes all whitespace characters from string s
void stripBlanks(char * s)
{
  char * w = s;
  while (*w != '\0')
  {
    while (*w == ' ') // erase blank
    {
      char * u = w;
      while (*u != '\0') // shift everything left one char
      {
        *u = *(u+1);
        u++;
      }
    }
    w++;
  }
}

int compareLoadList(const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int loadCommaList(const char * filename, int * numListEntries, int ** listEntries, int offset)
{
   // comma-separated text file of fixed vertices
  FILE * fin;
  fin = fopen(filename,"r");
  if (!fin)
  {
    printf("Error: could not open file %s.\n",filename);
    return 1;
  }

  *numListEntries = 0;

  char s[4096];
  while (fgets(s,4096,fin) != NULL)
  {
    stripBlanks(s);

    char * pch;
    pch = strtok (s,",");
    while ((pch != NULL) && (isdigit(*pch)))
    {
      (*numListEntries)++;
      pch = strtok (NULL, ",");
    }
  }

  *listEntries = (int*) malloc (sizeof(int) * (*numListEntries));

  rewind(fin);

  (*numListEntries) = 0;

  while (fgets(s,4096,fin) != NULL)
  {
    stripBlanks(s);
    char * pch;
    pch = strtok (s,",");
    while ((pch != NULL) && (isdigit(*pch)))
    {
      (*listEntries)[*numListEntries] = atoi(pch) - offset;
      (*numListEntries)++;
      pch = strtok (NULL, ",");
    }
  }

  // sort the list entries
  qsort ((*listEntries), *numListEntries, sizeof(int), compareLoadList);

  fclose(fin);

  return 0;
}

int loadListComparator(const void * a, const void * b)
{
  if (*(int*)a < *(int*)b)
    return -1;

  if (*(int*)a == *(int*)b)
    return 0;

  return 1;
}
