#include <QImageReader>
#include <QImageWriter>
#include <set>
#include "simOBJ.h"
#include "string_formatter.h"
#include "print_macro.h"

void SplitGroup(char *string, int *tri) {
  char buffer[1024];
  char* p;
  char* n;
  int slash = 0;
  int i = 0;
  int c = 0;

  strcpy(buffer, string);

  p = n = buffer;

  while (*p != 0) {
    if (*p == '/')
      slash++;
    p++;
  }

  if (slash < 2) {
    tri[0] = -1;
    tri[1] = -1;
    tri[2] = -1;
    return;
  }

  p = buffer;

  while (*n != 0) {
    while (*n != '/' && *n != 0) {
      c++;
      n++;
    }
    if (*n != 0) {
      *n = 0;
      n++;
    }

    if (c != 0)
      tri[i++] = atoi(p);
    else
      tri[i++] = -1;

    p = n;
    c = 0;
  }
}


bool isNumber(char character) {
  if (((int)character >= (int)'0' && (int)character <= (int)'9') || character == '.')
    return true;
  else
    return false;
}


bool isNumber(char *string) {
  int len = strlen(string);
  for (int i = 0; i < len; i++) {
    if (isNumber(string[i]))
      continue;
    else
      return false;
  }
  return true;
}


void chop(char *buffer) {
  int len = strlen(buffer);
  int i = len - 1;

  while (buffer[i] == '\n') {
    buffer[i] = 0;
    i--;
  }
}


Texture::Texture(char *filename) {
  QImageReader img_reader;
  QString texture_file(filename);
  img_reader.setFileName(texture_file);
  strcpy(this->filename, filename);
  if (img_reader.canRead()) {
    QImage img = img_reader.read().convertToFormat(QImage::Format_RGB888);
    ASSERT(img.isNull() == false);
    //        img = img.scaled(512, 512);
    img = img.mirrored(false, true);
    glGenTextures(1, &textureName);
    glBindTexture(GL_TEXTURE_2D, textureName);
    static int i = 0;
    //    if (i == 10)
    //    std::fill_n(img.bits(), img.width() * img.height() * 3, (i * 40) % 256);
    QImageWriter writer(QString("/Users/dj/%1.bmp").arg(i), "bmp");
    writer.setQuality(100);
    writer.write(img);
    ++i;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
                 img.width(), img.height(),
                 0, GL_RGB, GL_UNSIGNED_BYTE, img.bits());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    this->width = img.width();
    this->height = img.height();
    P(img.width(), img.height(), textureName);
  } else {
    auto msg = dj::Format("texture file %z can not be read.", filename);
    L(msg);
  }
}


MTL_File::MTL_File(const char *filename, const char* path) {
  FILE* fp;
  int read_state = DEFAULT;
  (void) read_state;
  int nTokens;
  int mat = -1;
  char** tokens;
  char buffer[1024];

  nMaterials = 0;

  //First pass
  fp = fopen(filename, "r");

  if (fp == NULL)
    return;

  while (fgets(buffer, 1024, fp))	{
    chop(buffer);
    tokens = SplitString(buffer, &nTokens);

    if (tokens == NULL)
      continue;

    if (strcmp(tokens[0], "newmtl") == 0)
      nMaterials++;

    delete tokens;
  }

  fclose(fp);

  materials = new Material[nMaterials];

  for (int i = 0; i < nMaterials; i++) {
    strcpy(materials[i].map_Ka, "");
    strcpy(materials[i].map_Ks, "");
    strcpy(materials[i].map_Kd, "");
    strcpy(materials[i].map_Ns, "");
    strcpy(materials[i].map_d, "");
    strcpy(materials[i].map_bump, "");
    materials[i].t_Ka = NULL;
    materials[i].t_Ks = NULL;
    materials[i].t_Kd = NULL;
    materials[i].t_Ns = NULL;
    materials[i].t_d = NULL;
    materials[i].t_bump = NULL;
    materials[i].hasTexture = false;
  }

  fp = fopen(filename, "r");

  if (fp == NULL)
    return;

  while (fgets(buffer, 1024, fp)) {
    chop(buffer);
    tokens = SplitString(buffer, &nTokens);

    if (tokens == NULL)
      continue;

    if (strcmp(tokens[0], "newmtl") == 0 && nTokens >= 2) {
      //printf("newmtl found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      strcpy(materials[++mat].name, tokens[1]);
    } else if (strcmp(tokens[0], "Ka") == 0 && nTokens >= 4) {
      //printf("Ka found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      materials[mat].Ka[0] = atof(tokens[1]);
      materials[mat].Ka[1] = atof(tokens[2]);
      materials[mat].Ka[2] = atof(tokens[3]);
    } else if (strcmp(tokens[0], "Kd") == 0 && nTokens >= 4) {
      //printf("Kd found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      materials[mat].Kd[0] = atof(tokens[1]);
      materials[mat].Kd[1] = atof(tokens[2]);
      materials[mat].Kd[2] = atof(tokens[3]);
    } else if (strcmp(tokens[0], "Ks") == 0 && nTokens >= 4) {
      //printf("Ks found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      materials[mat].Ks[0] = atof(tokens[1]);
      materials[mat].Ks[1] = atof(tokens[2]);
      materials[mat].Ks[2] = atof(tokens[3]);
    } else if ((strcmp(tokens[0], "Tr") == 0 || strcmp(tokens[0], "d") == 0) && nTokens >= 2) {
      //printf("Tr found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      materials[mat].Tr = atof(tokens[1]);
    } else if (strcmp(tokens[0], "Ns") == 0 && nTokens >= 2) {
      //printf("Ns found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      materials[mat].Ns = atof(tokens[1]);
    } else if (strcmp(tokens[0], "illum") == 0 && nTokens >= 2) {
      //printf("illum found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      materials[mat].illum = atof(tokens[1]);
    } else if (strcmp(tokens[0], "map_Kd") == 0 && nTokens >= 2) {
      //printf("map_Kd found\n");
      //for(int i=0; i<nTokens; i++)
      //        printf("\tToken: %s\n", tokens[1]);
      std::string texture_pic = dj::Format("%z/%z", path, tokens[1]);
      P(texture_pic);
      strcpy(materials[mat].map_Kd, texture_pic.c_str());
      materials[mat].t_Kd = new Texture(materials[mat].map_Kd);
      materials[mat].hasTexture = true;
    } else if (strcmp(tokens[0], "map_Ks") == 0 && nTokens >= 2) {
      //printf("map_Kd found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      std::string texture_pic = dj::Format("%z/%z", path, tokens[1]);
      P(texture_pic);
      strcpy(materials[mat].map_Ks, texture_pic.c_str());
      materials[mat].t_Ks = new Texture(materials[mat].map_Ks);
      materials[mat].hasTexture = true;
    } else if (strcmp(tokens[0], "map_Ka") == 0 && nTokens >= 2) {
      //printf("map_Kd found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      std::string texture_pic = dj::Format("%z/%z", path, tokens[1]);
      P(texture_pic);
      strcpy(materials[mat].map_Ka, texture_pic.c_str());
      materials[mat].t_Ka = new Texture(materials[mat].map_Ka);
      materials[mat].hasTexture = true;
    } else if ((strcmp(tokens[0], "map_d") == 0 || strcmp(tokens[0], "map_Tr") == 0) && nTokens >= 2) {
      //printf("map_Kd found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      std::string texture_pic = dj::Format("%z/%z", path, tokens[1]);
      P(texture_pic);
      strcpy(materials[mat].map_d, texture_pic.c_str());
      materials[mat].t_d = new Texture(materials[mat].map_d);
      materials[mat].hasTexture = true;
    } else if (strcmp(tokens[0], "map_Ns") == 0 && nTokens >= 2) {
      //printf("map_Kd found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      std::string texture_pic = dj::Format("%z/%z", path, tokens[1]);
      P(texture_pic);
      strcpy(materials[mat].map_Ns, texture_pic.c_str());
      materials[mat].t_Ns = new Texture(materials[mat].map_Ns);
      materials[mat].hasTexture = true;
    } else if ((strcmp(tokens[0], "map_bump") == 0 || strcmp(tokens[0], "bump") == 0) && nTokens >= 2) {
      //printf("map_Kd found\n");
      //for(int i=0; i<nTokens; i++)
      //	printf("\tToken: %s\n", tokens[i]);
      std::string texture_pic = dj::Format("%z/%z", path, tokens[1]);
      P(texture_pic);
      strcpy(materials[mat].map_bump, texture_pic.c_str());
      materials[mat].t_bump = new Texture(materials[mat].map_bump);
      materials[mat].hasTexture = true;
    } else {
      //printf("Unknown command ignored: %s\n", tokens[0]);
    }
    delete tokens;
  }
  fclose(fp);
}


void MTL_File::PrintDetails() {
  for (int i = 0; i < nMaterials; i++) {
    printf("Name:\t%s\n", materials[i].name);
    printf("\tAmbient:   %f %f %f\n", materials[i].Ka[0], materials[i].Ka[1], materials[i].Ka[2]);
    printf("\tDiffuse:   %f %f %f\n", materials[i].Kd[0], materials[i].Kd[1], materials[i].Kd[2]);
    printf("\tSpecular:  %f %f %f\n", materials[i].Ks[0], materials[i].Ks[1], materials[i].Ks[2]);

    printf("\tShininess: %f\n", materials[i].Ns);
    printf("\tAlpha:     %f\n", materials[i].Tr);
    printf("\tTexture (ambient):   %s\n", materials[i].map_Ka);
    printf("\tTexture (diffuse):   %s\n", materials[i].map_Kd);
    printf("\tTexture (specular):  %s\n", materials[i].map_Ks);
    printf("\tTexture (shininess): %s\n", materials[i].map_Ns);
    printf("\tTexture (alpha):     %s\n", materials[i].map_d);
    printf("\tTexture (bump map):  %s\n", materials[i].map_bump);
  }
}


Material *MTL_File::FindMaterialByName(char *name) {
  for (int i = 0; i < nMaterials; i++) {
    if (strcmp(name, materials[i].name) == 0)
      return &materials[i];
  }
  return NULL;
}


OBJ_Model::OBJ_Model(const char *filename, const char* path) {
  char buffer[1024];
  char** tokens;
  int nTokens;
  int loopState1 = DEFAULT;
  int loopState2 = DEFAULT;
  int loopState3 = DEFAULT;
  (void) loopState3;
  FILE* fp;
  Material* currentMat = NULL;

  int vertex = 0;
  int triangle = 0;
  int normal = 0;
  int tex = 0;
  int group = 0;
  int matf = 0;
  int object = 0;
  int gtriangle = 0;
  int otriangle = 0;

  strcpy(objectFileName, filename);

  nVertices = nTriangles = nNormals = nTexCoords = nObjects = nGroups = nMaterialLibs = 0;

  //Count number of groups.
  fp = fopen(filename, "r");

  if (fp == NULL) {
    puts("File doesn't exist...");
    return;
  }

  while (fgets(buffer, 1024, fp)) {
    chop(buffer);
    tokens = SplitString(buffer, &nTokens);

    if (tokens == NULL)
      continue;

    if (strcmp(tokens[0], "g") == 0) {
    } else if (strcmp(tokens[0], "o") == 0) {
      nObjects++;
    } else if (strcmp(tokens[0], "mtllib") == 0) {
      nMaterialLibs++;
    } else if (strcmp(tokens[0], "usemtl") == 0) {
      nGroups++;
    }
    delete tokens;
  }
  fclose(fp);

  groupArray = (OBJ_Group*)malloc(sizeof(OBJ_Group) * nGroups);
  objectArray = (OBJ_Object*)malloc(sizeof(OBJ_Object) * nObjects);
  mtlArray = (MTL_File**)malloc(sizeof(MTL_File*)*nMaterialLibs);

  for (int i = 0; i < nGroups; i++)
    groupArray[i].nTriangles = 0;
  for (int i = 0; i < nObjects; i++)
    objectArray[i].nTriangles = 0;

  //Count number of vertices, triangles, normals and texture coordinates.
  fp = fopen(filename, "r");

  if (fp == NULL) {
    puts("File doesn't exist...");
    return;
  }

  while (fgets(buffer, 1024, fp)) {
    chop(buffer);
    tokens = SplitString(buffer, &nTokens);

    if (tokens == NULL)
      continue;

    if (strcmp(tokens[0], "v") == 0) {
      //Vertex command.
      loopState1 = DEFAULT;
      loopState2 = DEFAULT;
      loopState3 = DEFAULT;
      nVertices++;
    } else if (strcmp(tokens[0], "vt") == 0) {
      //Vertex texture command.
      loopState1 = DEFAULT;
      loopState2 = DEFAULT;
      loopState3 = DEFAULT;
      nTexCoords++;
    } else if (strcmp(tokens[0], "vn") == 0) {
      //Vertex normal command.
      loopState1 = DEFAULT;
      loopState2 = DEFAULT;
      loopState3 = DEFAULT;
      nNormals++;
    } else if (strcmp(tokens[0], "f") == 0) {
      //Face command.
      if (loopState1 == GROUP)
        groupArray[group - 1].nTriangles++;
      if (loopState2 == OBJECT)
        objectArray[object - 1].nTriangles++;
      nTriangles++;
    } else if (strcmp(tokens[0], "g") == 0) {
      //Group command.
      loopState1 = GROUP;
      group++;
    } else if (strcmp(tokens[0], "o") == 0) {
      //Object command.
      loopState2 = OBJECT;
      object++;
    } else {
      //puts("Unknown command ignored");
    }
    delete tokens;
  }
  fclose(fp);
  fp = NULL;

  //Allocate arrays
  vertexArray		= new OBJ_Vertex[nVertices * 2];
  normalArray		= new OBJ_Normal[nNormals * 2];
  texCoordArray	= new OBJ_TexCoord[nTexCoords * 2];
  triangleArray	= new OBJ_Triangle[nTriangles * 2];

  for (int i = 0; i < nTriangles; i++) {
    triangleArray[i].Vertex[0]		= -1;
    triangleArray[i].Vertex[1]		= -1;
    triangleArray[i].Vertex[2]		= -1;
    triangleArray[i].TexCoord[0]	= -1;
    triangleArray[i].TexCoord[1]	= -1;
    triangleArray[i].TexCoord[1]	= -1;
    triangleArray[i].Normal[0]		= -1;
    triangleArray[i].Normal[1]		= -1;
    triangleArray[i].Normal[2]		= -1;
    triangleArray[i].material		= NULL;
  }

  for (int i = 0; i < nGroups; i++) {
    groupArray[i].triangleArray = new OBJ_Triangle[groupArray[i].nTriangles * 2];
    for (int j = 0; j < groupArray[i].nTriangles; j++) {
      groupArray[i].triangleArray[j].Vertex[0]	= -1;
      groupArray[i].triangleArray[j].Vertex[1]	= -1;
      groupArray[i].triangleArray[j].Vertex[2]	= -1;
      groupArray[i].triangleArray[j].TexCoord[0]	= -1;
      groupArray[i].triangleArray[j].TexCoord[1]	= -1;
      groupArray[i].triangleArray[j].TexCoord[2]	= -1;
      groupArray[i].triangleArray[j].Normal[0]	= -1;
      groupArray[i].triangleArray[j].Normal[1]	= -1;
      groupArray[i].triangleArray[j].Normal[2]	= -1;
      groupArray[i].triangleArray[j].material		= NULL;
    }
  }

  for (int i = 0; i < nObjects; i++) {
    objectArray[i].triangleArray = new OBJ_Triangle[objectArray[i].nTriangles * 2];
    for (int j = 0; j < objectArray[i].nTriangles; j++) {
      objectArray[i].triangleArray[j].Vertex[0]	= -1;
      objectArray[i].triangleArray[j].Vertex[1]	= -1;
      objectArray[i].triangleArray[j].Vertex[2]	= -1;
      objectArray[i].triangleArray[j].TexCoord[0]	= -1;
      objectArray[i].triangleArray[j].TexCoord[1]	= -1;
      objectArray[i].triangleArray[j].TexCoord[2]	= -1;
      objectArray[i].triangleArray[j].Normal[0]	= -1;
      objectArray[i].triangleArray[j].Normal[1]	= -1;
      objectArray[i].triangleArray[j].Normal[2]	= -1;
      objectArray[i].triangleArray[j].material	= NULL;
    }
  }

  group = 0;
  object = 0;
  currentMat = NULL;
  loopState1 = DEFAULT;
  loopState2 = DEFAULT;
  loopState3 = DEFAULT;

  //Now actually read the information
  fp = fopen(filename, "r");

  if (fp == NULL) {
    puts("File doesn't exist...");
    return;
  }

  std::set<int> ids;
  while (fgets(buffer, 1024, fp)) {
    chop(buffer);
    tokens = SplitString(buffer, &nTokens);

    if (tokens == NULL)
      continue;

    if (strcmp(tokens[0], "v") == 0) {
      //Vertex command.
      loopState1 = DEFAULT;
      loopState2 = DEFAULT;
      loopState3 = DEFAULT;
      if (nTokens >= 4) {
        vertexArray[vertex].X = atof(tokens[1]);
        vertexArray[vertex].Y = atof(tokens[2]);
        vertexArray[vertex++].Z = atof(tokens[3]);
      }
    } else if (strcmp(tokens[0], "vt") == 0) {
      //Vertex texture command.
      loopState1 = DEFAULT;
      loopState2 = DEFAULT;
      loopState3 = DEFAULT;
      if (nTokens >= 3) {
        texCoordArray[tex].U = atof(tokens[1]);
        texCoordArray[tex++].V = atof(tokens[2]);
      }
    } else if (strcmp(tokens[0], "vn") == 0) {
      //Vertex normal command.
      loopState1 = DEFAULT;
      loopState2 = DEFAULT;
      loopState3 = DEFAULT;
      if (nTokens >= 4) {
        normalArray[normal].X = atof(tokens[1]);
        normalArray[normal].Y = atof(tokens[2]);
        normalArray[normal++].Z = atof(tokens[3]);
      }
    } else if (strcmp(tokens[0], "f") == 0) {
      //Face command.
      //Fix this!
      if (nTokens >= 4) {
        if (loopState1 == GROUP) {
          sscanf(tokens[1], "%d/%d/%d",
                 &(groupArray[group - 1].triangleArray[gtriangle].Vertex[0]),
                 &(groupArray[group - 1].triangleArray[gtriangle].TexCoord[0]),
                 &(groupArray[group - 1].triangleArray[gtriangle].Normal[0]));
          sscanf(tokens[2], "%d/%d/%d",
                 &(groupArray[group - 1].triangleArray[gtriangle].Vertex[1]),
                 &(groupArray[group - 1].triangleArray[gtriangle].TexCoord[1]),
                 &(groupArray[group - 1].triangleArray[gtriangle].Normal[1]));
          sscanf(tokens[3], "%d/%d/%d",
                 &(groupArray[group - 1].triangleArray[gtriangle].Vertex[2]),
                 &(groupArray[group - 1].triangleArray[gtriangle].TexCoord[2]),
                 &(groupArray[group - 1].triangleArray[gtriangle].Normal[2]));
          groupArray[group - 1].triangleArray[gtriangle].material = currentMat;
          gtriangle++;
        }
        if (loopState2 == OBJECT) {
          sscanf(tokens[1], "%d/%d/%d",
                 &(objectArray[object - 1].triangleArray[otriangle].Vertex[0]),
                 &(objectArray[object - 1].triangleArray[otriangle].TexCoord[0]),
                 &(objectArray[object - 1].triangleArray[otriangle].Normal[0]));
          sscanf(tokens[2], "%d/%d/%d",
                 &(objectArray[object - 1].triangleArray[otriangle].Vertex[1]),
                 &(objectArray[object - 1].triangleArray[otriangle].TexCoord[1]),
                 &(objectArray[object - 1].triangleArray[otriangle].Normal[1]));
          sscanf(tokens[3], "%d/%d/%d",
                 &(objectArray[object - 1].triangleArray[otriangle].Vertex[2]),
                 &(objectArray[object - 1].triangleArray[otriangle].TexCoord[2]),
                 &(objectArray[object - 1].triangleArray[otriangle].Normal[2]));
          objectArray[object - 1].triangleArray[otriangle].material = currentMat;
          otriangle++;
        }
        sscanf(tokens[1], "%d/%d/%d",
               &(triangleArray[triangle].Vertex[0]),
               &(triangleArray[triangle].TexCoord[0]),
               &(triangleArray[triangle].Normal[0]));
        sscanf(tokens[2], "%d/%d/%d",
               &(triangleArray[triangle].Vertex[1]),
               &(triangleArray[triangle].TexCoord[1]),
               &(triangleArray[triangle].Normal[1]));
        sscanf(tokens[3], "%d/%d/%d",
               &(triangleArray[triangle].Vertex[2]),
               &(triangleArray[triangle].TexCoord[2]),
               &(triangleArray[triangle].Normal[2]));
        triangleArray[triangle].material = currentMat;
        if (currentMat->hasTexture && currentMat->t_Kd) {
          ids.insert(currentMat->t_Kd->textureName);
        }
        triangle++;
      }
    } else if (strcmp(tokens[0], "o") == 0) {
      //Object command.
      otriangle = 0;
      if (nTokens >= 2) {
        loopState2 = OBJECT;
        strcpy(objectArray[object].name, tokens[1]);
        objectArray[object++].name[strlen(tokens[1])] = '\0';
      }
    } else if (strcmp(tokens[0], "g") == 0) {
      //Group command.
      //      gtriangle = 0;
      //      if (nTokens >= 2) {
      //        loopState1 = GROUP;
      //        strcpy(groupArray[group].name, tokens[1]);
      //        groupArray[group++].name[strlen(tokens[1])] = '\0';
      //      }
    } else if (strcmp(tokens[0], "usemtl") == 0) {
      //Material command.
      Material* mat;
      for (int i = 0; i < nMaterialLibs; i++) {
        if ((mat = mtlArray[i]->FindMaterialByName(tokens[1])) != NULL)
          currentMat = mat;
        else {
          currentMat = NULL;
          auto msg = dj::Format("material \"%z\"not found", tokens[1]);
          ASSERT(false, L(msg));
        }
      }

      gtriangle = 0;
      if (nTokens >= 2) {
        loopState1 = GROUP;
        auto str = dj::Format("%d", group);
        strcpy(groupArray[group].name, str.c_str());
        groupArray[group].mtl = currentMat;//name[strlen(tokens[1])] = '\0';
        groupArray[group++].name[strlen(tokens[1])] = '\0';
      }
    } else if (strcmp(tokens[0], "mtllib") == 0) {
      //Material lib command.
      std::string mtl_file = dj::Format("%z/%z", path, tokens[1]);
      P(mtl_file);
      mtlArray[matf++] = new MTL_File(mtl_file.c_str(), path);
    } else {
      //puts("Unknown command ignored");
    }
    delete tokens;
  }
  fclose(fp);
  fp = NULL;

  P(ids.size(), nTriangles);
  for (auto& str : ids) {
    P(str);
  }
  {
    std::set<int> id;
    int cnt = 0;
    for (int i = 0; i < this->nTriangles; i++) {
      if (this->triangleArray[i].material) {
        if (this->triangleArray[i].material->hasTexture) {
          cnt++;
          id.insert(this->triangleArray[i].material->t_Kd->textureName);
          //          glBindTexture(GL_TEXTURE_2D, this->triangleArray->material->t_Kd->textureName);
        }
      }
    }
    P(id.size(), cnt);
    for (int i : id) {
      P(i);
    }
    //    exit(0);
  }
  //  exit(0);
  if (nNormals <= 0)
    this->RecalculateNormals();

  this->CalculateBoundingBoxes();
}


void OBJ_Model::CalculateBoundingBoxes() {
  for (int i = 0; i < nGroups; i++) {
    float minX, minY, minZ, maxX, maxY, maxZ;
    float x, y, z;
    minX = 100000.0;
    minY = 100000.0;
    minZ = 100000.0;
    maxX = -100000.0;
    maxY = -100000.0;
    maxZ = -100000.0;
    for (int j = 0; j < groupArray[i].nTriangles; j++) {
      for (int k = 0; k < 3; k++) {
        x = vertexArray[groupArray[i].triangleArray[j].Vertex[k] - 1].X;
        y = vertexArray[groupArray[i].triangleArray[j].Vertex[k] - 1].Y;
        z = vertexArray[groupArray[i].triangleArray[j].Vertex[k] - 1].Z;

        if (x < minX)	minX = x;
        if (y < minY)	minY = y;
        if (z < minZ)	minZ = z;
        if (x > maxX)	maxX = x;
        if (y > maxY)	maxY = y;
        if (z > maxZ)	maxZ = z;
      }
      groupArray[i].boundingBox.minX = minX;
      groupArray[i].boundingBox.minY = minY;
      groupArray[i].boundingBox.minZ = minZ;
      groupArray[i].boundingBox.maxX = maxX;
      groupArray[i].boundingBox.maxY = maxY;
      groupArray[i].boundingBox.maxZ = maxZ;
    }
  }

  for (int i = 0; i < nObjects; i++) {
    float minX, minY, minZ, maxX, maxY, maxZ;
    float x, y, z;
    minX = 100000.0;
    minY = 100000.0;
    minZ = 100000.0;
    maxX = -100000.0;
    maxY = -100000.0;
    maxZ = -100000.0;
    for (int j = 0; j < objectArray[i].nTriangles; j++) {
      for (int k = 0; k < 3; k++) {
        x = vertexArray[objectArray[i].triangleArray[j].Vertex[k] - 1].X;
        y = vertexArray[objectArray[i].triangleArray[j].Vertex[k] - 1].Y;
        z = vertexArray[objectArray[i].triangleArray[j].Vertex[k] - 1].Z;

        if (x < minX)	minX = x;
        if (y < minY)	minY = y;
        if (z < minZ)	minZ = z;
        if (x > maxX)	maxX = x;
        if (y > maxY)	maxY = y;
        if (z > maxZ)	maxZ = z;
      }
      objectArray[i].boundingBox.minX = minX;
      objectArray[i].boundingBox.minY = minY;
      objectArray[i].boundingBox.minZ = minZ;
      objectArray[i].boundingBox.maxX = maxX;
      objectArray[i].boundingBox.maxY = maxY;
      objectArray[i].boundingBox.maxZ = maxZ;
    }
  }
}


OBJ_Normal *OBJ_Model::CaculateNormal(OBJ_Triangle * tri) {
  OBJ_Normal* normal = new OBJ_Normal();

  //Calculate vectors
  float var1_x = vertexArray[tri->Vertex[1] - 1].X - vertexArray[tri->Vertex[0] - 1].X;
  float var1_y = vertexArray[tri->Vertex[1] - 1].Y - vertexArray[tri->Vertex[0] - 1].Y;
  float var1_z = vertexArray[tri->Vertex[1] - 1].Z - vertexArray[tri->Vertex[0] - 1].Z;

  float var2_x = vertexArray[tri->Vertex[2] - 1].X - vertexArray[tri->Vertex[0] - 1].X;
  float var2_y = vertexArray[tri->Vertex[2] - 1].Y - vertexArray[tri->Vertex[0] - 1].Y;
  float var2_z = vertexArray[tri->Vertex[2] - 1].Z - vertexArray[tri->Vertex[0] - 1].Z;

  //Get cross product of vectors
  normal->X = (var1_y * var2_z) - (var2_y * var1_z);
  normal->Y = (var1_z * var2_x) - (var2_z * var1_x);
  normal->Z = (var1_x * var2_y) - (var2_x * var1_y);

  //Normalize final vector
  float vLen = sqrt((normal->X * normal->X) + (normal->Y * normal->Y) + (normal->Z * normal->Z));

  normal->X = normal->X / vLen;
  normal->Y = normal->Y / vLen;
  normal->Z = normal->Z / vLen;

  return normal;
}


void OBJ_Model::RecalculateNormals() {
  OBJ_Normal* normal;
  nNormals = nTriangles * 3;
  normalArray = new OBJ_Normal[nNormals];

  for (int i = 0; i < nTriangles; i++) {
    normal = this->CaculateNormal(&this->triangleArray[i]);
    this->normalArray[i].X = normal->X;
    this->normalArray[i].Y = normal->Y;
    this->normalArray[i].Z = normal->Z;
    this->triangleArray[i].Normal[0] = i + 1;
    this->triangleArray[i].Normal[1] = i + 1;
    this->triangleArray[i].Normal[2] = i + 1;
    this->normalArray[this->triangleArray[i].Normal[0]].X += normal->X;
    this->normalArray[this->triangleArray[i].Normal[0]].X /= 2;
    this->normalArray[this->triangleArray[i].Normal[0]].Y += normal->Y;
    this->normalArray[this->triangleArray[i].Normal[0]].Y /= 2;
    this->normalArray[this->triangleArray[i].Normal[0]].Z += normal->Z;
    this->normalArray[this->triangleArray[i].Normal[0]].Z /= 2;
    this->normalArray[this->triangleArray[i].Normal[1]].X += normal->X;
    this->normalArray[this->triangleArray[i].Normal[1]].X /= 2;
    this->normalArray[this->triangleArray[i].Normal[1]].Y += normal->Y;
    this->normalArray[this->triangleArray[i].Normal[1]].Y /= 2;
    this->normalArray[this->triangleArray[i].Normal[1]].Z += normal->Z;
    this->normalArray[this->triangleArray[i].Normal[1]].Z /= 2;
    this->normalArray[this->triangleArray[i].Normal[2]].X += normal->X;
    this->normalArray[this->triangleArray[i].Normal[2]].X /= 2;
    this->normalArray[this->triangleArray[i].Normal[2]].Y += normal->Y;
    this->normalArray[this->triangleArray[i].Normal[2]].Y /= 2;
    this->normalArray[this->triangleArray[i].Normal[2]].Z += normal->Z;
    this->normalArray[this->triangleArray[i].Normal[2]].Z /= 2;

    delete normal;
  }
}


void OBJ_Model::PrintDetails(bool verbose) {
  PrintTextBox('*', '*', "Details for Model: %s", this->objectFileName);
  fprintf(stderr, "-Number of Vertices:  %8d\n", nVertices);
  fprintf(stderr, "-Number of Normals:   %8d\n", nNormals);
  fprintf(stderr, "-Number of TexCoords: %8d\n", nTexCoords);
  fprintf(stderr, "-Total Triangles:     %8d\n", nTriangles);
  PrintTextBox('-', '|', "Objects (%d)", nObjects);
  for (int i = 0; i < nObjects; i++) {
    fprintf(stderr, "\t+Name: %s\n", objectArray[i].name);
    fprintf(stderr, "\t\t+Triangles: %2d\n", objectArray[i].nTriangles);
    if (verbose) {
      for (int j = 0; j < objectArray[i].nTriangles; j++) {
        printf("\t\t\t+Triangle %8d\n", j);
        printf("\t\t\t\t-V1: %f %f %f\n", vertexArray[objectArray[i].triangleArray[j].Vertex[0] - 1].X, vertexArray[objectArray[i].triangleArray[j].Vertex[0] - 1].Y, vertexArray[objectArray[i].triangleArray[j].Vertex[0] - 1].Z);
        printf("\t\t\t\t-N1: %f %f %f\n", normalArray[objectArray[i].triangleArray[j].Normal[0] - 1].X, normalArray[objectArray[i].triangleArray[j].Normal[0] - 1].Y, normalArray[objectArray[i].triangleArray[j].Normal[0] - 1].Z);
        printf("\t\t\t\t-T1: %f %f\n", texCoordArray[objectArray[i].triangleArray[j].TexCoord[0] - 1].U, texCoordArray[objectArray[i].triangleArray[j].TexCoord[0] - 1].V);
        printf("\t\t\t\t-V2: %f %f %f\n", vertexArray[objectArray[i].triangleArray[j].Vertex[1] - 1].X, vertexArray[objectArray[i].triangleArray[j].Vertex[1] - 1].Y, vertexArray[objectArray[i].triangleArray[j].Vertex[1] - 1].Z);
        printf("\t\t\t\t-N2: %f %f %f\n", normalArray[objectArray[i].triangleArray[j].Normal[1] - 1].X, normalArray[objectArray[i].triangleArray[j].Normal[1] - 1].Y, normalArray[objectArray[i].triangleArray[j].Normal[1] - 1].Z);
        printf("\t\t\t\t-T2: %f %f\n", texCoordArray[objectArray[i].triangleArray[j].TexCoord[1] - 1].U, texCoordArray[objectArray[i].triangleArray[j].TexCoord[1] - 1].V);
        printf("\t\t\t\t-V3: %f %f %f\n", vertexArray[objectArray[i].triangleArray[j].Vertex[2] - 1].X, vertexArray[objectArray[i].triangleArray[j].Vertex[2] - 1].Y, vertexArray[objectArray[i].triangleArray[j].Vertex[2] - 1].Z);
        printf("\t\t\t\t-N3: %f %f %f\n", normalArray[objectArray[i].triangleArray[j].Normal[2] - 1].X, normalArray[objectArray[i].triangleArray[j].Normal[2] - 1].Y, normalArray[objectArray[i].triangleArray[j].Normal[2] - 1].Z);
        printf("\t\t\t\t-T3: %f %f\n", texCoordArray[objectArray[i].triangleArray[j].TexCoord[2] - 1].U, texCoordArray[objectArray[i].triangleArray[j].TexCoord[2] - 1].V);
        printf("\t\t\t\t-M:  %s\n", objectArray[i].triangleArray[j].material->name);
      }
    }
    fprintf(stderr, "\t\tBounding Box: \n");
    fprintf(stderr, "\t\t\tMin X: %f\n", objectArray[i].boundingBox.minX);
    fprintf(stderr, "\t\t\tMin Y: %f\n", objectArray[i].boundingBox.minY);
    fprintf(stderr, "\t\t\tMin Z: %f\n", objectArray[i].boundingBox.minZ);
    fprintf(stderr, "\t\t\tMax X: %f\n", objectArray[i].boundingBox.maxX);
    fprintf(stderr, "\t\t\tMax Y: %f\n", objectArray[i].boundingBox.maxY);
    fprintf(stderr, "\t\t\tMax Z: %f\n", objectArray[i].boundingBox.maxZ);
  }
  PrintTextBox('-', '|', "Groups (%d)", nGroups);
  for (int i = 0; i < nGroups; i++) {
    P(i);
    fprintf(stderr, "\t+Name: %s\n", groupArray[i].name);
    fprintf(stderr, "\t\t+Triangles: %2d\n", groupArray[i].nTriangles);
    if (verbose) {
      for (int j = 0; j < groupArray[i].nTriangles; j++) {
        printf("\t\t\t+Triangle %8d\n", j);
        printf("\t\t\t\t-V1: %f %f %f\n", vertexArray[groupArray[i].triangleArray[j].Vertex[0] - 1].X, vertexArray[groupArray[i].triangleArray[j].Vertex[0] - 1].Y, vertexArray[groupArray[i].triangleArray[j].Vertex[0] - 1].Z);
        printf("\t\t\t\t-N1: %f %f %f\n", normalArray[groupArray[i].triangleArray[j].Normal[0] - 1].X, normalArray[groupArray[i].triangleArray[j].Normal[0] - 1].Y, normalArray[groupArray[i].triangleArray[j].Normal[0] - 1].Z);
        printf("\t\t\t\t-T1: %f %f\n", texCoordArray[groupArray[i].triangleArray[j].TexCoord[0] - 1].U, texCoordArray[groupArray[i].triangleArray[j].TexCoord[0] - 1].V);
        printf("\t\t\t\t-V2: %f %f %f\n", vertexArray[groupArray[i].triangleArray[j].Vertex[1] - 1].X, vertexArray[groupArray[i].triangleArray[j].Vertex[1] - 1].Y, vertexArray[groupArray[i].triangleArray[j].Vertex[1] - 1].Z);
        printf("\t\t\t\t-N2: %f %f %f\n", normalArray[groupArray[i].triangleArray[j].Normal[1] - 1].X, normalArray[groupArray[i].triangleArray[j].Normal[1] - 1].Y, normalArray[groupArray[i].triangleArray[j].Normal[1] - 1].Z);
        printf("\t\t\t\t-T2: %f %f\n", texCoordArray[groupArray[i].triangleArray[j].TexCoord[1] - 1].U, texCoordArray[groupArray[i].triangleArray[j].TexCoord[1] - 1].V);
        printf("\t\t\t\t-V3: %f %f %f\n", vertexArray[groupArray[i].triangleArray[j].Vertex[2] - 1].X, vertexArray[groupArray[i].triangleArray[j].Vertex[2] - 1].Y, vertexArray[groupArray[i].triangleArray[j].Vertex[2] - 1].Z);
        printf("\t\t\t\t-N3: %f %f %f\n", normalArray[groupArray[i].triangleArray[j].Normal[2] - 1].X, normalArray[groupArray[i].triangleArray[j].Normal[2] - 1].Y, normalArray[groupArray[i].triangleArray[j].Normal[2] - 1].Z);
        printf("\t\t\t\t-T3: %f %f\n", texCoordArray[groupArray[i].triangleArray[j].TexCoord[2] - 1].U, texCoordArray[groupArray[i].triangleArray[j].TexCoord[2] - 1].V);
        printf("\t\t\t\t-M:  %s\n", groupArray[i].triangleArray[j].material->name);
      }
    }
    fprintf(stderr, "\t\tBounding Box: \n");
    fprintf(stderr, "\t\t\tMin X: %f\n", groupArray[i].boundingBox.minX);
    fprintf(stderr, "\t\t\tMin Y: %f\n", groupArray[i].boundingBox.minY);
    fprintf(stderr, "\t\t\tMin Z: %f\n", groupArray[i].boundingBox.minZ);
    fprintf(stderr, "\t\t\tMax X: %f\n", groupArray[i].boundingBox.maxX);
    fprintf(stderr, "\t\t\tMax Y: %f\n", groupArray[i].boundingBox.maxY);
    fprintf(stderr, "\t\t\tMax Z: %f\n", groupArray[i].boundingBox.maxZ);
  }
  KK;
  PrintTextBox('-', '|', "Mtllibs (%d)", nMaterialLibs);
  for (int i = 0; i < nMaterialLibs; i++) {
    mtlArray[i]->PrintDetails();
  }
  KK;
  printf("\n");
}


void OBJ_Model::Draw(int mode) {
  bool useInternalMaterials = mode & INTERNAL_MATERIALS;
  bool useInternalTextures = mode & INTERNAL_TEXTURES;

  //  std::set<int> id;
  static int i = 0;
  glBindTexture(GL_TEXTURE_2D, i % 13);
  ++i;
  for (int g = 0; g < nGroups; ++g) {
    P(g, nGroups);
    OBJ_Group& group = groupArray[g];
    if (group.mtl && group.mtl->hasTexture) {
      Material* mtl = group.mtl;
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, group.mtl->t_Ka->textureName);
      float Ka[4] = {mtl->Ka[0],
                     mtl->Ka[1],
                     mtl->Ka[2],
                     mtl->Tr
                    };

      float Ks[4] = {mtl->Ks[0],
                     mtl->Ks[1],
                     mtl->Ks[2],
                     mtl->Tr
                    };

      float Kd[4] = {mtl->Kd[0],
                     mtl->Kd[1],
                     mtl->Kd[2],
                     mtl->Tr
                    };
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ka);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Ks);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kd);
      glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, mtl->Ns);
    } else {
      glDisable(GL_TEXTURE_2D);
    }
    glBegin(GL_TRIANGLES);
    for (int t = 0; t < group.nTriangles; ++t) {
      OBJ_Triangle& tri = group.triangleArray[t];
      for (int i = 0; i < 3; ++i) {
        glTexCoord2f(this->texCoordArray[tri.TexCoord[i] - 1].U,
                     this->texCoordArray[tri.TexCoord[i] - 1].V);
        glNormal3f(this->normalArray[tri.Normal[i] - 1].X,
                   this->normalArray[tri.Normal[i] - 1].Y,
                   this->normalArray[tri.Normal[i] - 1].Z);
        glVertex3f(this->vertexArray[tri.Vertex[i] - 1].X,
                   this->vertexArray[tri.Vertex[i] - 1].Y,
                   this->vertexArray[tri.Vertex[i] - 1].Z);
      }

    }
    glEnd();
  }
  KK;
#if 0
  glBegin(GL_TRIANGLES);
  for (int i = 0; i < this->nTriangles; i++) {
    if (this->triangleArray[i].material) {
      if (useInternalMaterials) {
        float Ka[4] = {this->triangleArray[i].material->Ka[0],
                       this->triangleArray[i].material->Ka[1],
                       this->triangleArray[i].material->Ka[2],
                       this->triangleArray[i].material->Tr
                      };
        float Ks[4] = {this->triangleArray[i].material->Ks[0],
                       this->triangleArray[i].material->Ks[1],
                       this->triangleArray[i].material->Ks[2],
                       this->triangleArray[i].material->Tr
                      };
        float Kd[4] = {this->triangleArray[i].material->Kd[0],
                       this->triangleArray[i].material->Kd[1],
                       this->triangleArray[i].material->Kd[2],
                       this->triangleArray[i].material->Tr
                      };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ka);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Ks);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kd);
        glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, this->triangleArray[i].material->Ns);
      }
      if (useInternalTextures && this->triangleArray[i].material->hasTexture) {
        //        id.insert(this->triangleArray[i].material->t_Kd->textureName);
        glBindTexture(GL_TEXTURE_2D, this->triangleArray[i].material->t_Kd->textureName);
      }
    }
    if (useInternalTextures && this->triangleArray[i].material && this->triangleArray[i].material->hasTexture) {
      glTexCoord2f(this->texCoordArray[this->triangleArray[i].TexCoord[0] - 1].U,
                   this->texCoordArray[this->triangleArray[i].TexCoord[0] - 1].V);
    }
    glNormal3f(this->normalArray[this->triangleArray[i].Normal[0] - 1].X,
               this->normalArray[this->triangleArray[i].Normal[0] - 1].Y,
               this->normalArray[this->triangleArray[i].Normal[0] - 1].Z);
    glVertex3f(this->vertexArray[this->triangleArray[i].Vertex[0] - 1].X,
               this->vertexArray[this->triangleArray[i].Vertex[0] - 1].Y,
               this->vertexArray[this->triangleArray[i].Vertex[0] - 1].Z);
    if (useInternalTextures && this->triangleArray[i].material && this->triangleArray[i].material->hasTexture) {
      glTexCoord2f(this->texCoordArray[this->triangleArray[i].TexCoord[1] - 1].U,
                   this->texCoordArray[this->triangleArray[i].TexCoord[1] - 1].V);
    }
    glNormal3f(this->normalArray[this->triangleArray[i].Normal[1] - 1].X,
               this->normalArray[this->triangleArray[i].Normal[1] - 1].Y,
               this->normalArray[this->triangleArray[i].Normal[1] - 1].Z);
    glVertex3f(this->vertexArray[this->triangleArray[i].Vertex[1] - 1].X,
               this->vertexArray[this->triangleArray[i].Vertex[1] - 1].Y,
               this->vertexArray[this->triangleArray[i].Vertex[1] - 1].Z);
    if (useInternalTextures && this->triangleArray[i].material && this->triangleArray[i].material->hasTexture) {
      glTexCoord2f(this->texCoordArray[this->triangleArray[i].TexCoord[2] - 1].U,
                   this->texCoordArray[this->triangleArray[i].TexCoord[2] - 1].V);
    }
    glNormal3f(this->normalArray[this->triangleArray[i].Normal[2] - 1].X,
               this->normalArray[this->triangleArray[i].Normal[2] - 1].Y, this->normalArray[this->triangleArray[i].Normal[2] - 1].Z);
    glVertex3f(this->vertexArray[this->triangleArray[i].Vertex[2] - 1].X,
               this->vertexArray[this->triangleArray[i].Vertex[2] - 1].Y,
               this->vertexArray[this->triangleArray[i].Vertex[2] - 1].Z);
  }
  glEnd();
#endif
}


void OBJ_Model::DrawGroup(char *groupName, int mode) {
  OBJ_Group* group = NULL;
  bool useInternalMaterials = mode & INTERNAL_MATERIALS;
  bool useInternalTextures = mode & INTERNAL_TEXTURES;

  for (int i = 0; i < this->nGroups; i++) {
    if (strcmp(this->groupArray[i].name, groupName) == 0)	{
      group = groupArray + i;
      break;
    }
  }

  if (group == NULL)
    return;

  AABoundingBox* bb = &group->boundingBox;
  glTranslatef(-bb->maxX, -bb->maxY, -bb->maxZ);

  glBegin(GL_TRIANGLES );
  for (int i = 0; i < group->nTriangles; i++) {
    if (group->triangleArray[i].material) {
      if (useInternalMaterials) {
        float Ka[4] = {group->triangleArray[i].material->Ka[0],
                       group->triangleArray[i].material->Ka[1],
                       group->triangleArray[i].material->Ka[2],
                       group->triangleArray[i].material->Tr
                      };
        float Ks[4] = {group->triangleArray[i].material->Ks[0],
                       group->triangleArray[i].material->Ks[1],
                       group->triangleArray[i].material->Ks[2],
                       group->triangleArray[i].material->Tr
                      };
        float Kd[4] = {group->triangleArray[i].material->Kd[0],
                       group->triangleArray[i].material->Kd[1],
                       group->triangleArray[i].material->Kd[2],
                       group->triangleArray[i].material->Tr
                      };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ka);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Ks);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kd);
        glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, group->triangleArray[i].material->Ns);
      }
      if (useInternalTextures && group->triangleArray[i].material->hasTexture)
        glBindTexture(GL_TEXTURE_2D, group->triangleArray->material->t_Kd->textureName);
    }
    if (useInternalTextures && group->triangleArray[i].material && group->triangleArray[i].material->hasTexture)
      glTexCoord2f(this->texCoordArray[group->triangleArray[i].TexCoord[0] - 1].U,
                   this->texCoordArray[group->triangleArray[i].TexCoord[0] - 1].V);
    glNormal3f(this->normalArray[group->triangleArray[i].Normal[0] - 1].X,
               this->normalArray[group->triangleArray[i].Normal[0] - 1].Y,
               this->normalArray[group->triangleArray[i].Normal[0] - 1].Z);
    glVertex3f(this->vertexArray[group->triangleArray[i].Vertex[0] - 1].X,
               this->vertexArray[group->triangleArray[i].Vertex[0] - 1].Y,
               this->vertexArray[group->triangleArray[i].Vertex[0] - 1].Z);
    if (useInternalTextures && group->triangleArray[i].material && group->triangleArray[i].material->hasTexture)
      glTexCoord2f(this->texCoordArray[group->triangleArray[i].TexCoord[1] - 1].U,
                   this->texCoordArray[group->triangleArray[i].TexCoord[1] - 1].V);
    glNormal3f(this->normalArray[group->triangleArray[i].Normal[1] - 1].X,
               this->normalArray[group->triangleArray[i].Normal[1] - 1].Y,
               this->normalArray[group->triangleArray[i].Normal[1] - 1].Z);
    glVertex3f(this->vertexArray[group->triangleArray[i].Vertex[1] - 1].X,
               this->vertexArray[group->triangleArray[i].Vertex[1] - 1].Y,
               this->vertexArray[group->triangleArray[i].Vertex[1] - 1].Z);
    if (useInternalTextures && group->triangleArray[i].material && group->triangleArray[i].material->hasTexture)
      glTexCoord2f(this->texCoordArray[group->triangleArray[i].TexCoord[2] - 1].U,
                   this->texCoordArray[group->triangleArray[i].TexCoord[2] - 1].V);
    glNormal3f(this->normalArray[group->triangleArray[i].Normal[2] - 1].X,
               this->normalArray[group->triangleArray[i].Normal[2] - 1].Y,
               this->normalArray[group->triangleArray[i].Normal[2] - 1].Z);
    glVertex3f(this->vertexArray[group->triangleArray[i].Vertex[2] - 1].X,
               this->vertexArray[group->triangleArray[i].Vertex[2] - 1].Y,
               this->vertexArray[group->triangleArray[i].Vertex[2] - 1].Z);
  }
  glEnd();
}


void OBJ_Model::DrawObject(char *objectName, int mode) {
  OBJ_Object* object = NULL;
  bool useInternalMaterials = mode & INTERNAL_MATERIALS;
  bool useInternalTextures = mode & INTERNAL_TEXTURES;

  for (int i = 0; i < this->nObjects; i++) {
    if (strcmp(this->objectArray[i].name, objectName) == 0)	{
      object = groupArray + i;
      break;
    }
  }

  if (object == NULL)
    return;

  AABoundingBox* bb = &object->boundingBox;
  glTranslatef(-bb->maxX, -bb->maxY, -bb->maxZ);

  glBegin(GL_TRIANGLES );
  for (int i = 0; i < object->nTriangles; i++) {
    if (object->triangleArray[i].material) {
      if (useInternalMaterials) {
        float Ka[4] = {object->triangleArray[i].material->Ka[0], object->triangleArray[i].material->Ka[1], object->triangleArray[i].material->Ka[2], object->triangleArray[i].material->Tr};
        float Ks[4] = {object->triangleArray[i].material->Ks[0], object->triangleArray[i].material->Ks[1], object->triangleArray[i].material->Ks[2], object->triangleArray[i].material->Tr};
        float Kd[4] = {object->triangleArray[i].material->Kd[0], object->triangleArray[i].material->Kd[1], object->triangleArray[i].material->Kd[2], object->triangleArray[i].material->Tr};
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ka);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Ks);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kd);
        glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, object->triangleArray[i].material->Ns);
      }
      if (useInternalTextures && object->triangleArray[i].material->hasTexture)
        glBindTexture(GL_TEXTURE_2D, object->triangleArray->material->t_Kd->textureName);
    }
    if (useInternalTextures && object->triangleArray[i].material && object->triangleArray[i].material->hasTexture)
      glTexCoord2f(this->texCoordArray[object->triangleArray[i].TexCoord[0] - 1].U, this->texCoordArray[object->triangleArray[i].TexCoord[0] - 1].V);
    glNormal3f(this->normalArray[object->triangleArray[i].Normal[0] - 1].X, this->normalArray[object->triangleArray[i].Normal[0] - 1].Y, this->normalArray[object->triangleArray[i].Normal[0] - 1].Z);
    glVertex3f(this->vertexArray[object->triangleArray[i].Vertex[0] - 1].X, this->vertexArray[object->triangleArray[i].Vertex[0] - 1].Y, this->vertexArray[object->triangleArray[i].Vertex[0] - 1].Z);
    if (useInternalTextures && object->triangleArray[i].material && object->triangleArray[i].material->hasTexture)
      glTexCoord2f(this->texCoordArray[object->triangleArray[i].TexCoord[1] - 1].U, this->texCoordArray[object->triangleArray[i].TexCoord[1] - 1].V);
    glNormal3f(this->normalArray[object->triangleArray[i].Normal[1] - 1].X, this->normalArray[object->triangleArray[i].Normal[1] - 1].Y, this->normalArray[object->triangleArray[i].Normal[1] - 1].Z);
    glVertex3f(this->vertexArray[object->triangleArray[i].Vertex[1] - 1].X, this->vertexArray[object->triangleArray[i].Vertex[1] - 1].Y, this->vertexArray[object->triangleArray[i].Vertex[1] - 1].Z);
    if (useInternalTextures && object->triangleArray[i].material && object->triangleArray[i].material->hasTexture)
      glTexCoord2f(this->texCoordArray[object->triangleArray[i].TexCoord[2] - 1].U, this->texCoordArray[object->triangleArray[i].TexCoord[2] - 1].V);
    glNormal3f(this->normalArray[object->triangleArray[i].Normal[2] - 1].X, this->normalArray[object->triangleArray[i].Normal[2] - 1].Y, this->normalArray[object->triangleArray[i].Normal[2] - 1].Z);
    glVertex3f(this->vertexArray[object->triangleArray[i].Vertex[2] - 1].X, this->vertexArray[object->triangleArray[i].Vertex[2] - 1].Y, this->vertexArray[object->triangleArray[i].Vertex[2] - 1].Z);
  }
  glEnd();
}


AABoundingBox *OBJ_Model::getGroupBoundingBox(char *groupName) {
  for (int i = 0; i < nGroups; i++) {
    if (strcmp(groupArray[i].name, groupName) == 0)
      return &groupArray[i].boundingBox;
  }

  return NULL;
}


AABoundingBox *OBJ_Model::getObjectBoundingBox(char *objectName) {
  for (int i = 0; i < nObjects; i++) {
    if (strcmp(objectArray[i].name, objectName) == 0)
      return &objectArray[i].boundingBox;
  }

  return NULL;
}
