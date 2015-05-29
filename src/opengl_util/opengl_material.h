#ifndef OPENGLMATERIAL_H
#define OPENGLMATERIAL_H
#include "opengl_header.h"

class OpenGLMaterial {
  typedef float Real;
public:
  enum MaterialType {
    kEmerald,
    kJade,
    kObsidian,
    kPearl,
    kRuby,
    kTurquoise,
    kBrass,
    kBronze,
    kChrome,
    kCopper,
    kGold,
    kSilver,
    kBlackPlastic,
    kCyanPlastic,
    kGreenPlastic,
    kRedPlastic,
    kWhitePlastic,
    kYellowPlastic,
    kBlackRubber,
    kCyanRubber,
    kGreenRubber,
    kRedRubber,
    kWhiteRubber,
    kYellowRubber,
    kMaterialNum
  };

  OpenGLMaterial(const Real* ambient = nullptr,
                 const Real* diffuse = nullptr,
                 const Real* specular = nullptr,
                 Real shininess = Real(128));
  OpenGLMaterial(const OpenGLMaterial& material);
  const OpenGLMaterial& operator=(const OpenGLMaterial& material);
  ~OpenGLMaterial();
  void Use(GLenum material_face = GL_FRONT_AND_BACK) const;
  const Real* ambient() const { return ambient_; }
  const Real* diffuse() const { return diffuse_; }
  const Real* specular() const { return specular_; }
  Real shininess() const { return shininess_; }
  static OpenGLMaterial CreateMaterial(int material_type);
private:
  void Construct(const Real *ambient, const Real *diffuse, const Real *specular, Real shininess);
  Real ambient_[4];
  Real diffuse_[4];
  Real specular_[4];
  Real shininess_;
};

#endif // OPENGLMATERIAL_H
