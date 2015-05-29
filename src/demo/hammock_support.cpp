#include "hammock_support.h"
#include "reflective_object_renderer.h"
#include "opengl_helper.h"
#include "affine_transformer.h"
#include "triangle_mesh_creator.h"
#include "global.h"
#include "print_macro.h"
#include "string_formatter.h"
#include "binary_file_io.h"

HammockSupport::HammockSupport() {
  support_ = new ReflectiveObjectRenderer<double>(DATA_DIRECTORY "texture/rome_church", DATA_DIRECTORY "../shader/env_map");
  std::vector<double> expanded_verts;
  std::vector<double> expanded_normal;
  int triangle_num = 0;
//  if (0)
 // {
 //   const char* translate_commands[] = {
 //     "translate -0.98 0 -0.312",
 //     "translate -0.98 0 +0.312",
 //     "translate +0.98 0 +0.312",
 //     "translate +0.98 0 -0.312",
 //   };
    std::vector<int> tris;
    std::vector<double> verts, normals;
    //CreateCylinder(0.018, 1.0, 20, 20, verts, normals, tris);

    //Read verts binary
    std::string support_file_name = dj::Format("%z/verts_support.bin",GetDataFolder());
    BinaryFileReader in5(support_file_name.c_str());
    int numElems;
    in5.Read(&numElems, 1);
    verts.resize(numElems);
    in5.Read(&verts[0], numElems);
    P(numElems);
    //Read verts binary
    support_file_name = dj::Format("%z/normal_support.bin",GetDataFolder());
    BinaryFileReader in2(support_file_name.c_str());
    in2.Read(&numElems, 1);
    normals.resize(numElems);
    in2.Read(&normals[0], numElems);


    //Read tris binary
    support_file_name = dj::Format("%z/tris_support.bin",GetDataFolder());
    BinaryFileReader in3(support_file_name.c_str());
    in3.Read(&numElems, 1);
    tris.resize(numElems);
    in3.Read(&tris[0], numElems);

    triangle_num = numElems/3;

    /*std::vector<std::string> cmd;
    cmd.push_back("centerize 0 0 0");
    cmd.push_back("range 1.0");
    //cmd.push_back("centerize 0 0 0.1");
    AffineTransformer<double> transformer(cmd);
    transformer.Transform(&verts[0], int(verts.size() / 3));*/


    ExpandTriangleMesh(tris, verts, normals, expanded_verts, expanded_normal);




/*

    for (int i = 0; i < 4; ++i) {
      auto tmp_verts = verts;
      triangle_num += int(tris.size()) / 3;
      std::vector<std::string> cmd;
      cmd.push_back("max_y 0.51");
      cmd.push_back(translate_commands[i]);
      AffineTransformer<double> transformer(cmd);
      transformer.Transform(&tmp_verts[0], int(tmp_verts.size() / 3));

      std::vector<double> tmp_expanded_verts;
      std::vector<double> tmp_expanded_normal;
      ExpandTriangleMesh(tris, tmp_verts, normals, tmp_expanded_verts, tmp_expanded_normal);
      expanded_verts.insert(expanded_verts.end(), tmp_expanded_verts.begin(), tmp_expanded_verts.end());
      expanded_normal.insert(expanded_normal.end(), tmp_expanded_normal.begin(), tmp_expanded_normal.end());
    }
  }
//  if (0)
  {
    const char* translate_commands[] = {
      "translate -0.975 0.52 -0.312",
      "translate -0.98 0.52 +0.312",
      "translate +0.98 0.52 -0.312",
      "translate +0.98 0.52 +0.312",
    };
    std::vector<int> tris;
    std::vector<double> verts, normals;
    CreateCylinder(0.025, 0.07, 20, 20, verts, normals, tris);
//    CreateCylinder(0.025, 1.07, 20, 20, verts, normals, tris);
    {
      std::vector<std::string> cmd;
      cmd.push_back("rotate_x 90");
      AffineTransformer<double> transformer(cmd);
      transformer.Transform(&normals[0], int(normals.size() / 3));
    }
    for (int i = 0; i < 4; ++i) {
      auto tmp_verts = verts;
      triangle_num += int(tris.size()) / 3;
      std::vector<std::string> cmd;
      cmd.push_back("centerize");
      cmd.push_back("rotate_x 90");
      cmd.push_back(translate_commands[i]);
      AffineTransformer<double> transformer(cmd);
      transformer.Transform(&tmp_verts[0], int(tmp_verts.size() / 3));

      std::vector<double> tmp_expanded_verts;
      std::vector<double> tmp_expanded_normal;
      ExpandTriangleMesh(tris, tmp_verts, normals, tmp_expanded_verts, tmp_expanded_normal);
      expanded_verts.insert(expanded_verts.end(), tmp_expanded_verts.begin(), tmp_expanded_verts.end());
      expanded_normal.insert(expanded_normal.end(), tmp_expanded_normal.begin(), tmp_expanded_normal.end());
    }
  }
*/
  ASSERT(glGetError() == GL_NO_ERROR);
  support_->UpdatePosAndNormal(&expanded_verts[0], &expanded_normal[0], triangle_num);
  ASSERT(glGetError() == GL_NO_ERROR);
}

void HammockSupport::Render() {
  support_->Render();
}

HammockSupport::~HammockSupport() {
  delete support_;
}
