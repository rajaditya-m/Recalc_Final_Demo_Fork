#pragma once

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "cbmp.h"
#include "pconsole.h"
#include "strmanip.h"
#include "opengl_helper.h"

#define DEFAULT 0
#define GROUP	1
#define OBJECT	2
#define MTL		3

void SplitGroup(char* string, int* tri);

bool isNumber(char character);

bool isNumber(char* string);

void chop(char* buffer);

struct Texture {
  char filename[1024];
  int width, height;
//  CBmp* textureData;
  GLuint textureName;
  Texture(char* filename);
};

struct Material {
  float Ka[3];
  float Kd[3];
  float Ks[3];
  float Tr;
  float Ns;
  float illum;

  char map_Ka[1024];
  char map_Kd[1024];
  char map_Ks[1024];
  char map_d[1024];
  char map_Ns[1024];
  char map_bump[1024];

  Texture* t_Ka;
  Texture* t_Kd;
  Texture* t_Ks;
  Texture* t_d;
  Texture* t_Ns;
  Texture* t_bump;

  bool hasTexture;

  char name[1024];
};

class MTL_File {
private:
  int nMaterials;
  Material* materials;
public:
  MTL_File(const char* filename, const char *path);
  void PrintDetails();
  Material* FindMaterialByName(char* name);
};

struct AABoundingBox {
  float minX, minY, minZ;
  float maxX, maxY, maxZ;
};

struct OBJ_Vertex {
  float X, Y, Z;
};

typedef OBJ_Vertex OBJ_Normal;

struct OBJ_TexCoord {
  float U, V;
};

struct OBJ_Triangle {
  int Vertex[3];
  int Normal[3];
  int TexCoord[3];
  Material* material;
};

struct OBJ_Group {
  char name[1024];
  int nTriangles;
  Material* mtl;
  OBJ_Triangle* triangleArray;
  AABoundingBox boundingBox;
};

typedef struct OBJ_Group OBJ_Object;

#define INTERNAL_TEXTURES 1
#define INTERNAL_MATERIALS 2

class OBJ_Model {
private:
  int nVertices, nNormals, nTexCoords, nTriangles, nObjects, nGroups, nMaterials, nMaterialLibs;

  OBJ_Vertex		*vertexArray;
  OBJ_Normal		*normalArray;
  OBJ_TexCoord	*texCoordArray;

  OBJ_Triangle	*triangleArray;
  OBJ_Object		*objectArray;
  OBJ_Group		*groupArray;
  MTL_File		**mtlArray;

  char objectFileName[1024];

  OBJ_Normal* CaculateNormal(OBJ_Triangle *triangle);
  void CalculateBoundingBoxes();
public:
  OBJ_Model(const char* filename, const char *path);
  void PrintDetails(bool verbose = false);
  void RecalculateNormals();
  void Draw(int mode = 0);
  void DrawGroup(char* groupName, int mode = 0);
  void DrawObject(char* objectName, int mode = 0);
  AABoundingBox* getGroupBoundingBox(char* groupName);
  AABoundingBox* getObjectBoundingBox(char* objectName);
};

