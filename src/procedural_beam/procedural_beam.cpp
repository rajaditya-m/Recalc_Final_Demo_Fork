#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <iostream>
#include <set>
#include "procedural_beam.h"

struct TetraHedron {
  int a; int b; int c; int d;
  TetraHedron(int _a, int _b, int _c, int _d)
    : a(_a), b(_b), c(_c), d(_d)
  {};
};

#define ADD1(x) (x)+1

ProceduralBeam::ProceduralBeam(double lengthX, double lengthY, double lengthZ,
                               int numPointsAlongX, int numPointsAlongY, int numPointsAlongZ) {
  m_lengthX = lengthX;
  m_lengthY = lengthY;
  m_lengthZ = lengthZ;

  m_numPointsAlongX = numPointsAlongX;
  m_numPointsAlongY = numPointsAlongY;
  m_numPointsAlongZ = numPointsAlongZ;
}


ProceduralBeam::~ProceduralBeam(void) {
}

inline int ProceduralBeam::generateGlobalIndex(int r, int c, int z) {
  int numPointsPerCS = m_numPointsAlongX * m_numPointsAlongY;
  int res = (r * m_numPointsAlongY + c) + (z * numPointsPerCS);
  return res;
}

inline std::pair<int, int> ProceduralBeam::makeEdgeNicely(int x, int y) {
  if (x > y) {
    int t = x;
    x = y;
    y = t;
  }
  return std::pair<int, int>(x, y);
}

//Output
// NODE FILE - Done
// ELE FILE - Done
// VEG FILE - Done
// BOU FILE - Done
// SUR FILE - Done
// CROSS FILE - Done
// OBJ FILE - dONE
//EDGE FILE  - Done
void ProceduralBeam::generateTetsProcedurally(const char* suffix) {
  //compute the offsets here also
  double rowOffset = m_lengthX / (m_numPointsAlongX - 1);
  double colOffset = m_lengthY / (m_numPointsAlongY - 1);
  double lenOffset = m_lengthZ / (m_numPointsAlongZ - 1);

  //Some statistical information will be useful later
  int totPoints = m_numPointsAlongY * m_numPointsAlongX * m_numPointsAlongZ;

  //Our point data vector
  std::vector<Eigen::Vector3d> pointData;
  std::vector<TetraHedron> tets;
  std::string fSuffix(suffix);

  //Every cross section has 100 points
  double curLen = 0.0;
  for (int z = 0; z < m_numPointsAlongZ; z++) {
    double rowVal = 0.0;
    for (int r = 0; r < m_numPointsAlongX; r++) {
      double colVal = 0.0;
      for (int c = 0; c < m_numPointsAlongY; c++) {
        pointData.push_back(Eigen::Vector3d(rowVal, colVal, curLen));
        colVal += colOffset;
      }
      rowVal += rowOffset;
    }
    curLen += lenOffset;
  }
  //  curLen -= lenOffset;
  //  for (int z = 0; z < m_numCrossSection; z++) {
  //    double rowVal = 0.0;
  //    for (int r = 0; r < m_pointsPerRow; r++) {
  //      double colVal = 0.0;
  //      for (int c = 0; c < m_pointsPerCol; c++) {
  //        pointData.push_back(Eigen::Vector3d(rowVal, colVal, curLen));
  //        colVal += colOffset;
  //      }
  //      rowVal += rowOffset;
  //    }
  //    curLen += lenOffset;
  //  }

  //Now generate the tetrahedrons procedurally
  int a, b, c, d, e, f;
  //  int tetCounter = 0;
  for (int z = 0; z < m_numPointsAlongZ - 1; z++) {
    for (int x = 0; x < m_numPointsAlongX - 1; x++) {
      for (int y = 0; y < m_numPointsAlongY - 1; y++) {
        int i1 = generateGlobalIndex(x, y, z);
        int i2 = generateGlobalIndex(x, y + 1, z);
        int i3 = generateGlobalIndex(x + 1, y + 1, z);
        int i4 = generateGlobalIndex(x + 1, y, z);

        int j1 = generateGlobalIndex(x, y, z + 1);
        int j2 = generateGlobalIndex(x, y + 1, z + 1);
        int j3 = generateGlobalIndex(x + 1, y + 1, z + 1);
        int j4 = generateGlobalIndex(x + 1, y, z + 1);

        a = i1; b = i2; c = i3;
        d = j1; e = j2; f = j3;

        tets.push_back(TetraHedron(a, b, c, d));
        tets.push_back(TetraHedron(b, c, e, d));
        tets.push_back(TetraHedron(e, c, f, d));

        //a = i1; b = i3; c = i4;
        //d = j1; e = j3; f = j4;
        //a = i3; b = i4; c = i1;
        //d = j3; e = j4; f = j1;
        a = i1; b = i4; c = i3;
        d = j1; e = j4; f = j3;

        tets.push_back(TetraHedron(a, b, c, d));
        //tets.push_back(TetraHedron(a,b,c,f));
        tets.push_back(TetraHedron(b, c, e, d));
        //tets.push_back(TetraHedron(b,f,c,d));
        tets.push_back(TetraHedron(e, c, f, d));
        //tets.push_back(TetraHedron(b,c,f,d));

      }
    }
  }
  //  for (int z = m_numCrossSection; z < (m_numCrossSection * 2) - 1; z++) {
  //    for (int r = 0; r < m_pointsPerRow - 1; r++) {
  //      for (int ci = 0; ci < m_pointsPerCol - 1; ci++) {
  //        int i1 = generateGlobalIndex(r, ci, z);
  //        int i2 = generateGlobalIndex(r, ci + 1, z);
  //        int i3 = generateGlobalIndex(r + 1, ci + 1, z);
  //        int i4 = generateGlobalIndex(r + 1, ci, z);

  //        int j1 = generateGlobalIndex(r, ci, z + 1);
  //        int j2 = generateGlobalIndex(r, ci + 1, z + 1);
  //        int j3 = generateGlobalIndex(r + 1, ci + 1, z + 1);
  //        int j4 = generateGlobalIndex(r + 1, ci, z + 1);

  //        a = i1; b = i2; c = i3;
  //        d = j1; e = j2; f = j3;

  //        tets.push_back(TetraHedron(a, b, c, d));
  //        tets.push_back(TetraHedron(b, c, e, d));
  //        tets.push_back(TetraHedron(e, c, f, d));

  //        //a = i1; b = i3; c = i4;
  //        //d = j1; e = j3; f = j4;
  //        //a = i3; b = i4; c = i1;
  //        //d = j3; e = j4; f = j1;
  //        a = i1; b = i4; c = i3;
  //        d = j1; e = j4; f = j3;

  //        tets.push_back(TetraHedron(a, b, c, d));
  //        //tets.push_back(TetraHedron(a,b,c,f));
  //        tets.push_back(TetraHedron(b, c, e, d));
  //        //tets.push_back(TetraHedron(b,f,c,d));
  //        tets.push_back(TetraHedron(e, c, f, d));
  //        //tets.push_back(TetraHedron(b,c,f,d));

  //      }
  //    }
  //  }

  //Generate the NODE FILE
  std::string nodeFileName = fSuffix + std::string(".node");
  std::ofstream nodeFile(nodeFileName.c_str());
  //  nodeFile << totPoints * 2 << " 3 0 0\n";
  nodeFile << totPoints << " 3 0 0\n";
  for (int i = 0; i < totPoints; i++) {
    nodeFile << i << " " << pointData[i].x() << " " << pointData[i].y() << " " << pointData[i].z() << "\n";
  }
  nodeFile.close();

  //Generate the ELE FILE Procedurally
  std::string eleFileName = fSuffix + std::string(".ele");
  std::ofstream eleFile(eleFileName.c_str());
  eleFile << tets.size() << " 4 0\n";
  for (int t = 0; t < int(tets.size()); t++) {
    eleFile << t << " " << tets[t].a << " " << tets[t].b << " " << tets[t].c << " " << tets[t].d << "\n";
  }
  eleFile.close();
#if 0
  //Generate the VEG FILE now
  std::string vegFileName = fSuffix + std::string(".veg");
  std::ofstream vegFile(vegFileName.c_str());
  vegFile << "*VERTICES\n";
  vegFile << "*INCLUDE " << nodeFileName << "\n\n";
  vegFile << "*ELEMENTS\n";
  vegFile << "TETS\n";
  vegFile << "*INCLUDE " << eleFileName << "\n\n";
  vegFile << "*MATERIAL defaultMat\n";
  vegFile << "ENU, 1.1E3, 1E6, 0.49\n\n";
  vegFile << "*REGION\n";
  vegFile << "allElements, defaultMat\n";
  vegFile.close();

  //Generate the OBJ File
  std::string objFileName = fSuffix + std::string(".obj");
  std::ofstream objFile(objFileName);
  for (int i = 0; i < totPoints * 2; i++) {
    objFile << "v " << pointData[i].x() << " " << pointData[i].y() << " " << pointData[i].z() << "\n";
  }
  objFile << "\n";
  for (int r = 0; r < m_numPointsAlongX - 1; r++) {
    for (int c = 0; c < m_numPointsAlongY - 1; c++) {
      objFile << "f " << ADD1(generateGlobalIndex(r, c, 0)) << " " << ADD1(generateGlobalIndex(r, c + 1, 0)) << " " << ADD1(generateGlobalIndex(r + 1, c + 1, 0)) << "\n";
      objFile << "f " << ADD1(generateGlobalIndex(r, c, 0)) << " " << ADD1(generateGlobalIndex(r + 1, c + 1, 0)) << " " << ADD1(generateGlobalIndex(r + 1, c, 0)) << "\n";

      objFile << "f " << ADD1(generateGlobalIndex(r, c, 2 * m_numPointsAlongZ - 1)) << " " << ADD1(generateGlobalIndex(r, c + 1, 2 * m_numPointsAlongZ - 1)) << " " << ADD1(generateGlobalIndex(r + 1, c + 1, 2 * m_numPointsAlongZ - 1)) << "\n";
      objFile << "f " << ADD1(generateGlobalIndex(r, c, 2 * m_numPointsAlongZ - 1)) << " " << ADD1(generateGlobalIndex(r + 1, c + 1, 2 * m_numPointsAlongZ - 1)) << " " << ADD1(generateGlobalIndex(r + 1, c, 2 * m_numPointsAlongZ - 1)) << "\n";
    }
  }
  for (int z = 0; z < (m_numPointsAlongZ * 2) - 1; z++) {
    for (int r = 0; r < m_numPointsAlongX - 1; r++) {
      objFile << "f " << ADD1(generateGlobalIndex(r, 0, z)) << " " << ADD1(generateGlobalIndex(r, 0, z + 1)) << " " << ADD1(generateGlobalIndex(r + 1, 0, z + 1)) << "\n";
      objFile << "f " << ADD1(generateGlobalIndex(r + 1, 0, z + 1)) << " " << ADD1(generateGlobalIndex(r + 1, 0, z)) << " " << ADD1(generateGlobalIndex(r, 0, z)) << "\n";

      objFile << "f " << ADD1(generateGlobalIndex(r, m_numPointsAlongY - 1, z)) << " " << ADD1(generateGlobalIndex(r, m_numPointsAlongY - 1, z + 1)) << " " << ADD1(generateGlobalIndex(r + 1, m_numPointsAlongY - 1, z + 1)) << "\n";
      objFile << "f " << ADD1(generateGlobalIndex(r + 1, m_numPointsAlongY - 1, z + 1)) << " " << ADD1(generateGlobalIndex(r + 1, m_numPointsAlongY - 1, z)) << " " << ADD1(generateGlobalIndex(r, m_numPointsAlongY - 1, z)) << "\n";
    }

    for (int c = 0; c < m_numPointsAlongY - 1; c++) {
      objFile << "f " << ADD1(generateGlobalIndex(0, c, z)) << " " << ADD1(generateGlobalIndex(0, c, z + 1)) << " " << ADD1(generateGlobalIndex(0, c + 1, z + 1)) << "\n";
      objFile << "f " << ADD1(generateGlobalIndex(0, c + 1, z + 1)) << " " << ADD1(generateGlobalIndex(0, c + 1, z)) << " " << ADD1(generateGlobalIndex(0, c, z)) << "\n";

      objFile << "f " << ADD1(generateGlobalIndex(m_numPointsAlongX - 1, c, z)) << " " << ADD1(generateGlobalIndex(m_numPointsAlongX - 1, c, z + 1)) << " " << ADD1(generateGlobalIndex(m_numPointsAlongX - 1, c + 1, z + 1)) << "\n";
      objFile << "f " << ADD1(generateGlobalIndex(m_numPointsAlongX - 1, c + 1, z + 1)) << " " << ADD1(generateGlobalIndex(m_numPointsAlongX - 1, c + 1, z)) << " " << ADD1(generateGlobalIndex(m_numPointsAlongX - 1, c, z)) << "\n";
    }
  }
  objFile.close();

  //Bou File Generation
  std::string bouFileName = fSuffix + std::string(".bou");
  std::ofstream bouFile(bouFileName);
  for (int z = 0; z < (m_numPointsAlongZ * 2); z++) {
    for (int r = 0; r < m_numPointsAlongX; r++) {
      bouFile << ADD1(generateGlobalIndex(r, 0, z)) << ",";
    }
    bouFile << "\n";
  }
  bouFile.close();

  //Surface File Generation
  std::string surFileName = fSuffix + std::string(".sur");
  std::ofstream surFile(surFileName);
  for (int z = 9; z < 15; z++) {
    for (int c = 1; c < m_numPointsAlongY; c++) {
      surFile << ADD1(generateGlobalIndex(0, c, z)) << "\n";
      surFile << ADD1(generateGlobalIndex(m_numPointsAlongX - 1, c, z)) << "\n";
    }
    for (int r = 1; r < m_numPointsAlongX - 1; r++) {
      surFile << ADD1(generateGlobalIndex(r, m_numPointsAlongY - 1, z)) << "\n";
    }
  }
  surFile.close();

  //Cross File Generation
  std::string crossFileName = fSuffix + std::string(".cross");
  std::ofstream crossFile(crossFileName);
  for (int z = 0; z < m_numPointsAlongZ; z++) {
    for (int r = 0; r < m_numPointsAlongX; r++) {
      for (int c = 0; c < m_numPointsAlongY; c++) {
        crossFile << ADD1(generateGlobalIndex(r, c, z)) << " " ;
      }
    }
    crossFile << "\n";
  }
  crossFile.close();

  //Edge file generation
  std::string edgeFileName = fSuffix + std::string(".edge");
  std::ofstream edgeFile(edgeFileName.c_str());
  std::set<std::pair<int, int> > edgeSet;
  for (int z = 0; z < m_numPointsAlongZ - 1; z++) {
    for (int r = 0; r < m_numPointsAlongX - 1; r++) {
      for (int ci = 0; ci < m_numPointsAlongY - 1; ci++) {
        int i1 = generateGlobalIndex(r, ci, z);
        int i2 = generateGlobalIndex(r, ci + 1, z);
        int i3 = generateGlobalIndex(r + 1, ci + 1, z);
        int i4 = generateGlobalIndex(r + 1, ci, z);

        int j1 = generateGlobalIndex(r, ci, z + 1);
        int j2 = generateGlobalIndex(r, ci + 1, z + 1);
        int j3 = generateGlobalIndex(r + 1, ci + 1, z + 1);
        int j4 = generateGlobalIndex(r + 1, ci, z + 1);

        //There will be a total of 24 edges

        //i1 and i2
        edgeSet.insert(makeEdgeNicely(i1, i2));
        //i2 and i3
        edgeSet.insert(makeEdgeNicely(i2, i3));
        //i3 and i4
        edgeSet.insert(makeEdgeNicely(i3, i4));
        //i4 and i1
        edgeSet.insert(makeEdgeNicely(i4, i1));
        //j1 and j2
        edgeSet.insert(makeEdgeNicely(j1, j2));
        //j2 and j3
        edgeSet.insert(makeEdgeNicely(j2, j3));
        //j3 and j4
        edgeSet.insert(makeEdgeNicely(j3, j4));
        //j4 and j1
        edgeSet.insert(makeEdgeNicely(j4, j1));
        // i1 and j1
        edgeSet.insert(makeEdgeNicely(i1, j1));
        // i2 and j2
        edgeSet.insert(makeEdgeNicely(i2, j2));
        // i3 and j3
        edgeSet.insert(makeEdgeNicely(i3, j3));
        // i4 and j4
        edgeSet.insert(makeEdgeNicely(i4, j4));

        //i1 and i3
        edgeSet.insert(makeEdgeNicely(i1, i3));
        //i2 and i4
        edgeSet.insert(makeEdgeNicely(i2, i4));
        //j1 and j3
        edgeSet.insert(makeEdgeNicely(j1, j3));
        //j2 and j4
        edgeSet.insert(makeEdgeNicely(j2, j4));
        //i3 and j2
        edgeSet.insert(makeEdgeNicely(i3, j2));
        //i2 and j3
        edgeSet.insert(makeEdgeNicely(i2, j3));
        //i4 and j1
        edgeSet.insert(makeEdgeNicely(i4, j1));
        //i1 and j4
        edgeSet.insert(makeEdgeNicely(i1, j4));
        //i1 and j2
        edgeSet.insert(makeEdgeNicely(i1, j2));
        //i2 and j1
        edgeSet.insert(makeEdgeNicely(i2, j1));
        //i4 and j3
        edgeSet.insert(makeEdgeNicely(i4, j3));
        //i3 and j4
        edgeSet.insert(makeEdgeNicely(i3, j4));

      }
    }
  }
  std::cout << "[INFO] Number of isotropic edges generated:" << edgeSet.size() << "\n";
  for (auto it = edgeSet.begin(); it != edgeSet.end(); it++) {
    edgeFile << ADD1(it->first) << " " << ADD1(it->second) << "\n";
  }
  edgeFile.close();
#endif
}

