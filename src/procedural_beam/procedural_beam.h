#pragma once
#include <utility>

class ProceduralBeam {
public:
  ProceduralBeam(double lengthX, double lengthY, double lengthZ,
                 int numPointsAlongX,
                 int numPointsAlongY,
                 int numPointsAlongZ);
  ~ProceduralBeam(void);

  void generateTetsProcedurally(const char* fname);

  //other helper functions
  inline int generateGlobalIndex(int r, int c, int z);
  inline std::pair<int, int> makeEdgeNicely(int x, int y);


private :
  double m_lengthX;
  double m_lengthY;
  double m_lengthZ;

  int m_numPointsAlongX;
  int m_numPointsAlongY;
  int m_numPointsAlongZ;
};

