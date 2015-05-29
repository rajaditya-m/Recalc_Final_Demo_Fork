#include "BLOCK_MATRIX_GRAPH.h"

namespace solver {

template <class T>
BLOCK_MATRIX_GRAPH<T>::BLOCK_MATRIX_GRAPH(std::vector<std::vector<int> > &topology, std::vector<int> &block_size, int additional_entry) {
  max_number = 32;
  max_e_number = 128;
  Av  = new MATRIX<T>[max_number];
  Ae  = new MATRIX<T>[max_e_number];
  E = new int[max_e_number * 2];
  VE  = new int[max_number][64];
  ve_number = new int[max_number];

  number = (int) topology.size();
  assert(number < max_number);
  e_number = 0;
  for (int v = 0; v < number; ++v) {
    for (int i = v + 1; i < number; ++i) {
      if (topology[v][i]) {
        E[e_number * 2 + 0] = v;
        E[e_number * 2 + 1] = i;
        e_number++;
      }
    }
  }
  e_number0 = e_number;

  topology = std::vector<std::vector<int> >(number, std::vector<int>(number, INT_MIN));
  for (int e = 0; e < e_number; ++e) {
    int v0 = E[e * 2 + 0];
    int v1 = E[e * 2 + 1];
    topology[v0][v1] = topology[v1][v0] = e;
  }
  for (int v = 0; v < number; ++v) {
    topology[v][v] = v;
  }

  assert(e_number < max_e_number);
  assert(number == block_size.size());
  size_ = 0;
  for (int v = 0; v < number; ++v) {
    int size = block_size[v] + additional_entry;
    Av[v].Set(0, size, size);
    size_ += size;
  }
  for (int e = 0; e < e_number; ++e) {
    int v0 = E[e * 2 + 0];
    int v1 = E[e * 2 + 1];
    Ae[e].Set(0, block_size[v0] + additional_entry, block_size[v1] + additional_entry);
  }
  BuildVE();
  BuildEdgeIndex();
}


template <class T>
BLOCK_MATRIX_GRAPH<T>::BLOCK_MATRIX_GRAPH(int block_num, std::vector<int> edges, std::vector<int> block_size) {
  max_number = 32;
  max_e_number = 128;
  Av  = new MATRIX<T>[max_number];
  Ae  = new MATRIX<T>[max_e_number];
  E = new int[max_e_number * 2];
  VE  = new int[max_number][64];
  ve_number = new int[max_number];

  number = block_num;
  assert(number < max_number);
  e_number = int(edges.size()) / 2;
  memcpy(E, &edges[0], sizeof(int) * edges.size());
  e_number0 = e_number;

  assert(e_number < max_e_number);
  assert(number == block_size.size());
  size_ = 0;
  for (int v = 0; v < number; ++v) {
    int size = block_size[v];
    Av[v].Set(0, size, size);
    size_ += size;
  }
  for (int e = 0; e < e_number; ++e) {
    int v0 = E[e * 2 + 0];
    int v1 = E[e * 2 + 1];
    Ae[e].Set(0, block_size[v0], block_size[v1]);
  }
  BuildVE();
  BuildEdgeIndex();
}

template <class T>
int BLOCK_MATRIX_GRAPH<T>::Cholesky_Decomposition() {
  MATRIX<T> temp_A;

  for (int i = 0; i < number; i++) {
    //First, process the vertex matrix (the diagonal)

    for (int eid = 0; eid < ve_number[i]; eid++) {
      int e = VE[i][eid];
      int j = (E[e * 2 + 0] == i) ? E[e * 2 + 1] : E[e * 2 + 0];
      if (j < i) {
        MATRIX<T>::Transpose_Multiply(Ae[e], Ae[e], temp_A);
        MATRIX<T>::Subtract(Av[i], temp_A, Av[i]);
      }
    }
    if (Av[i].Cholesky_Decomposition() == false) {
      return i;
    }
    //Now get its inverse for fast processing later
    Av[i].Upper_Inverse(temp_A);
    Av[i].Get(temp_A);

    //Second, process the existing edge matrix (the off-diagonal)
    for (int eid = 0; eid < ve_number[i]; eid++) {
      int ei = VE[i][eid];
      int k = (E[ei * 2 + 0] == i) ? E[ei * 2 + 1] : E[ei * 2 + 0];
      if (k < i) {
        for (int ekd = 0; ekd < ve_number[k]; ekd++) {
          int ek = VE[k][ekd];
          int j = (E[ek * 2 + 0] == k) ? E[ek * 2 + 1] : E[ek * 2 + 0];

          if (j > i) { //Unprocessed
            int ej = INT_MAX, ejd = INT_MAX;
            for (ejd = 0; ejd < ve_number[i]; ejd++) {
              ej = VE[i][ejd];
              if (E[ej * 2 + 0] == j || E[ej * 2 + 1] == j)  break;
            }
            if (ejd == ve_number[i]) { //a new edge
              ej = e_number++;
              E[ej * 2 + 0] = i;
              E[ej * 2 + 1] = j;
              VE[i][ve_number[i]++] = ej;
              VE[j][ve_number[j]++] = ej;
              Ae[ej].Set(0, Av[i].ni, Av[j].ni);
            }

            //Now the edges are: ei [k, i]; ek [k, j]; ej [i, j]
            MATRIX<T>::Transpose_Multiply(Ae[ei], Ae[ek], temp_A);
            MATRIX<T>::Subtract(Ae[ej], temp_A, Ae[ej]);
          }
        }
      }
    }

    //Finally, finish off-diagonal matrices
    for (int eid = 0; eid < ve_number[i]; eid++) {
      int e = VE[i][eid];
      int j = INT_MIN;
      if (E[e * 2 + 0] == i)  j = E[e * 2 + 1];
      else if (E[e * 2 + 1] == i)  j = E[e * 2 + 0];
      if (j > i) {
        //Av[i].Upper_Inverse(At1);
        MATRIX<T>::Transpose_Multiply(Av[i], Ae[e], temp_A);
        Ae[e].Get(temp_A);
      }
    }
  }

  //printf("E: %d\n", e_number);
  //printf("E8: %d, %d\n", E[16], E[17]);
  //Ae[8].Print();
  //Av[4].Upper_Inverse(temp_A);
  //temp_A.Print();
  return kSuccess;
}

template <class T>
int BLOCK_MATRIX_GRAPH<T>::SolveWithDecomposedMatrix(T b[], T x[]) {
  //Set up the indices in the x array
  int *index = new int[number + 1];
  index[0] = 0;
  for (int i = 1; i <= number; i++)
    index[i] = index[i - 1] + Av[i - 1].ni;

  memcpy(x, b, sizeof(T)*index[number]);

  //Lower solve: forward
  T temp[1024];
  for (int i = 0; i < number; i++) {
    for (int eid = 0; eid < ve_number[i]; eid++) {
      int e = VE[i][eid];
      int j = INT_MIN;
      if (E[e * 2 + 0] == i)  j = E[e * 2 + 1];
      else if (E[e * 2 + 1] == i)  j = E[e * 2 + 0];

      if (j < i) {
        MATRIX<T>::Transpose_Multiply(Ae[e], &x[index[j]], temp);
        for (int l = 0; l < Av[i].ni; l++) x[index[i] + l] -= temp[l];
      }
    }
    MATRIX<T>::Transpose_Multiply(Av[i], &x[index[i]], temp);
    for (int l = 0; l < Av[i].ni; l++) x[index[i] + l] = temp[l];
  }

  //Upper solve: backward
  for (int i = number - 1; i >= 0; i--) {
    for (int eid = 0; eid < ve_number[i]; eid++) {
      int e = VE[i][eid];
      int j = -100000;
      if (E[e * 2 + 0] == i)  j = E[e * 2 + 1];
      else if (E[e * 2 + 1] == i)  j = E[e * 2 + 0];

      if (j > i) {
        MATRIX<T>::Multiply(Ae[e], &x[index[j]], temp);
        for (int l = 0; l < Av[i].ni; l++) x[index[i] + l] -= temp[l];
      }
    }
    MATRIX<T>::Multiply(Av[i], &x[index[i]], temp);
    for (int l = 0; l < Av[i].ni; l++) x[index[i] + l] = temp[l];
  }

  delete[] index;
  return kSuccess;
}

template <class T>
int BLOCK_MATRIX_GRAPH<T>::Solve(T b[], T x[]) {
  //Prepare the matrix by decomposition
  int code = Cholesky_Decomposition();
  if (code != kSuccess) {
    return code;
  }

  //Set up the indices in the x array
  int *index = new int[number + 1];
  index[0] = 0;
  for (int i = 1; i <= number; i++)
    index[i] = index[i - 1] + Av[i - 1].ni;

  memcpy(x, b, sizeof(T)*index[number]);

  //Lower solve: forward
  T temp[1024];
  for (int i = 0; i < number; i++) {
    for (int eid = 0; eid < ve_number[i]; eid++) {
      int e = VE[i][eid];
      int j = INT_MIN;
      if (E[e * 2 + 0] == i)  j = E[e * 2 + 1];
      else if (E[e * 2 + 1] == i)  j = E[e * 2 + 0];

      if (j < i) {
        MATRIX<T>::Transpose_Multiply(Ae[e], &x[index[j]], temp);
        for (int l = 0; l < Av[i].ni; l++) x[index[i] + l] -= temp[l];
      }
    }
    MATRIX<T>::Transpose_Multiply(Av[i], &x[index[i]], temp);
    for (int l = 0; l < Av[i].ni; l++) x[index[i] + l] = temp[l];
  }

  //Upper solve: backward
  for (int i = number - 1; i >= 0; i--) {
    for (int eid = 0; eid < ve_number[i]; eid++) {
      int e = VE[i][eid];
      int j = -100000;
      if (E[e * 2 + 0] == i)  j = E[e * 2 + 1];
      else if (E[e * 2 + 1] == i)  j = E[e * 2 + 0];

      if (j > i) {
        MATRIX<T>::Multiply(Ae[e], &x[index[j]], temp);
        for (int l = 0; l < Av[i].ni; l++) x[index[i] + l] -= temp[l];
      }
    }
    MATRIX<T>::Multiply(Av[i], &x[index[i]], temp);
    for (int l = 0; l < Av[i].ni; l++) x[index[i] + l] = temp[l];
  }

  delete[] index;
  return kSuccess;
}

template <class T>
void BLOCK_MATRIX_GRAPH<T>::Multiply(T x[], T r[]) {
  //Set up the indices in the x array
  int *index = new int[number];
  index[0] = 0;
  for (int i = 1; i < number; i++)
    index[i] = index[i - 1] + Av[i - 1].ni;

  T temp[1024];
  for (int i = 0; i < number; i++) {
    //Multiply the diagonal part
    MATRIX<T>::Multiply(Av[i], &x[index[i]], &r[index[i]]);

    //Multiply the off-diagonal part
    for (int eid = 0; eid < ve_number[i]; eid++) {
      int e = VE[i][eid];
      int j = INT_MIN;
      if (E[e * 2 + 0] == i)  j = E[e * 2 + 1];
      else if (E[e * 2 + 1] == i)  j = E[e * 2 + 0];

      if (j < i) {
        MATRIX<T>::Transpose_Multiply(Ae[e], &x[index[j]], temp);
        for (int l = 0; l < Av[i].ni; l++) r[index[i] + l] += temp[l];
      } else {
        MATRIX<T>::Multiply(Ae[e], &x[index[j]], temp);
        for (int l = 0; l < Av[i].ni; l++) r[index[i] + l] += temp[l];
      }
    }
  }
}

template <class T>
void BLOCK_MATRIX_GRAPH<T>::BuildEdgeIndex() {
  vv2edge_ = std::vector<int>(number * number, INT_MIN);
  for (int i = 0; i < number; ++i) {
    vv2edge_[i * number + i] = i;
  }
  for (int e = 0; e < e_number0; ++e) {
    int v0 = E[e * 2 + 0];
    int v1 = E[e * 2 + 1];
    vv2edge_[v0 * number + v1] = e;
    vv2edge_[v1 * number + v0] = e;
  }
}

template class BLOCK_MATRIX_GRAPH<double>;
}
