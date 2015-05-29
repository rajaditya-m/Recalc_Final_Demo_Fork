#include "main_window.h"
#include <QDebug>
#include <QFileDialog>
#include <QApplication>
#include "input_handler.h"
#include "rainbow_color.h"
#include "global.h"
#include "open_gl_qt.h"
#include "print_macro.h"
#include "opengl_helper.h"
#if 1
#include "objMeshRender.h"
#include "objMesh.h"
class ObjMtlRender : public InputHandler {
public:
  ObjMtlRender(const char* obj_file)
    : obj_(obj_file)
    , renderer_(&obj_) {
    //    int textureMode = OBJMESHRENDER_GL_USEANISOTROPICFILTERING | OBJMESHRENDER_GL_USEMIPMAP | OBJMESHRENDER_GL_REPLACE;
    int textureMode = OBJMESHRENDER_GL_USEANISOTROPICFILTERING | OBJMESHRENDER_GL_USEMIPMAP | OBJMESHRENDER_LIGHTINGMODULATIONBIT;
    renderer_.loadTextures(textureMode);
  }
  virtual int Render() {
    glEnable(GL_LIGHTING);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glEnable(GL_NORMALIZE);
    //    glBegin(GL_QUADS);
    //    glNormal3f(0, 0, 1);
    //    glVertex3f(0, 0, 0);
    //    glVertex3f(0, 1, 0);
    //    glVertex3f(1, 1, 0);
    //    glVertex3f(1, 0, 0);
    //    glEnd();
    const float s = 0.00015;
    //    const float s = 0.51015;
    glScalef(s, s, s);
    renderer_.render(OBJMESHRENDER_TRIANGLES, OBJMESHRENDER_MATERIAL | OBJMESHRENDER_TEXTURE | OBJMESHRENDER_SMOOTH);
    //    renderer_.render(OBJMESHRENDER_TRIANGLES, OBJMESHRENDER_TEXTURE | OBJMESHRENDER_SMOOTH);
    //    exit(0);
    if (0) {
      //      float diffuse_color[4] = { 1, 1, 1, 0.6f};
      //      float dark[4] = {0, 0, 0, 1};
      //      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse_color);
      //      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, dark);
      glPushMatrix();
      //      glTranslatef(1 / s, 0, 0);
      //      P(obj_.getNumGroups());
      for (int g = 0; g < obj_.getNumGroups(); ++g) {
        const ObjMesh::Group* group = obj_.getGroupHandle(g);
        auto materialHandle = obj_.getMaterialHandle(group->getMaterialIndex());

        Vec3d Ka = materialHandle->getKa();
        Vec3d Kd = materialHandle->getKd();
        Vec3d Ks = materialHandle->getKs();

        float shininess = materialHandle->getShininess();
        float alpha = materialHandle->getAlpha();
        float ambient[4] = { (float)Ka[0], (float)Ka[1], (float)Ka[2], alpha };
        float diffuse[4] = { (float)Kd[0], (float)Kd[1], (float)Kd[2], alpha };
        float specular[4] = { (float)Ks[0], (float)Ks[1], (float)Ks[2], alpha };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
        if (materialHandle->hasTextureFilename()) {
          ObjMeshRender::Texture* textureHandle = renderer_.getTextureHandle(group->getMaterialIndex());
          glBindTexture(GL_TEXTURE_2D, textureHandle->getTexture());
          glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
          glEnable(GL_TEXTURE_2D);
        } else {
          glDisable(GL_TEXTURE_2D);
        }
        glBegin(GL_TRIANGLES);
        for (int f = 0; f < group->getNumFaces(); ++f) {
          const ObjMesh::Face* face = group->getFaceHandle(f);
          for (int i = 0; i < face->getNumVertices(); ++i) {
            const ObjMesh::Vertex* vert = face->getVertexHandle(i);
            Vec3d pos = obj_.getPosition(*vert);
            Vec3d normal = obj_.getNormal(*vert);
            Vec3d vtexture = obj_.getTextureCoordinate(*vert);
            glTexCoord2d(vtexture[0], vtexture[1]);
            //            dj::Normalize3(&normal[0]);
            //            normal *= 1 / s * 0.01;
            //            normal += pos;
            //          normal *= 0.1;
            Normal(&normal[0]);
            Vertex3v(&pos[0]);
            //          DrawArrow(&pos[0], &normal[0], true, 1.00, 0.01, 0.05);
            //          DrawArrow(&pos[0], &normal[0], true);
          }
        }
        glEnd();
        glDisable(GL_TEXTURE_2D);
        //      break;
      }
      //      exit(0);
      glPopMatrix();
    }
    glPopMatrix();
    //    glDisable(GL_BLEND);
    //    glDisable(GL_TEXTURE_2D);
    return kNotHandled;
  }

  ObjMesh obj_;
  ObjMeshRender renderer_;
};

void Init() {
  //  static ObjMtlRender render("/Users/dj/Dropbox/3D Model of Octopus/OCTOPUS/OCTOPUS/octopus.3.obj");
  //  static ObjMtlRender render("/Users/dj/Dropbox/3D Model of Octopus/OCTOPUS/OCTOPUS/OCTOPUS_.obj");
  static ObjMtlRender render("/Users/dj/Dropbox/3D Model of Octopus/OCTOPUS/OCTOPUS/octopus.3.obj");
  //    static ObjMtlRender render(DATA_DIRECTORY "obj/cube.obj");
  //    static ObjMtlRender render("/Users/dj/Downloads/capsule.obj");
  global::gl->InstallHandler(&render);
  //    glBindTexture(GL_TEXTURE_2D, 6);
}

int main(int argc, char * argv[]) {
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
#ifdef WINDOW_TITLE
#if defined(_WIN32) || defined(_WIN64)
  QString title(STRINGIZE_TOKEN(WINDOW_TITLE));
#else
  QString title(WINDOW_TITLE);
#endif
  w.setWindowTitle(title);
#endif
  w.show();
  w.hide();
  Init(); // Make sure OpenGL environment is setup before init
  w.show();
  return a.exec();
}
#else
#include "simOBJ.h"
class ObjMtlRender : public InputHandler {
public:
  ObjMtlRender(const char* obj_file, const char* path)
    : obj_(obj_file, path) {
    obj_.PrintDetails();
  }
  virtual int Render() {
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glPushMatrix();
    glColor3fv(kRed());
    const float s = 0.00015;
    glScalef(s, s, s);
    obj_.Draw(INTERNAL_MATERIALS | INTERNAL_TEXTURES);
    glPopMatrix();
    glDisable(GL_BLEND);
    glDisable(GL_TEXTURE_2D);
    return kNotHandled;
  }

  OBJ_Model obj_;
};

void Init() {
  //   static ObjMtlRender render("/Users/dj/Dropbox/3D Model of Octopus/OCTOPUS/OCTOPUS/octopus.3.obj",
  //                             "/Users/dj/Dropbox/3D Model of Octopus/OCTOPUS/OCTOPUS");
  static ObjMtlRender render("/Users/dj/3D Model of Octopus/OCTOPUS/OCTOPUS/octopus.3.obj",
                             "/Users/dj/3D Model of Octopus/OCTOPUS/OCTOPUS");
  //  static ObjMtlRender render(DATA_DIRECTORY "obj/cube.obj", DATA_DIRECTORY "obj");
  //  static ObjMtlRender render("/Users/dj/Downloads/capsule.obj", "/Users/dj/Downloads");
  global::gl->InstallHandler(&render);
  //    glBindTexture(GL_TEXTURE_2D, 6);
}

int main(int argc, char * argv[]) {
  QApplication a(argc, argv);
  setlocale(LC_NUMERIC, "C");
  MainWindow w;
  w.show();
  Init();
  return a.exec();
}
#endif


