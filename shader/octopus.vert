varying vec3 positionEye;
varying vec3 frag_normal;
varying vec3 LightDir;
varying vec3 ViewDir;

//attribute vec3 tangent;
attribute vec2 tex;
attribute vec3 position;
attribute vec3 normal;
//varying vec2 texCoord ;
//uniform sampler2D texture;
//uniform sampler2D bumpMap;
//uniform float test;

void main() {
  //  texCoord = tex;
//  positionEye = gl_ModelViewMatrix * gl_Vertex;
  positionEye = (gl_ModelViewMatrix * vec4(position, 1.0)).xyz;
  vec3 transformedNormal = gl_NormalMatrix * normal;
  vec3 normalEye = normalize(transformedNormal);
  frag_normal = normalEye;
//  gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
  gl_Position = gl_ModelViewProjectionMatrix * vec4(position, 1.0);
  gl_TexCoord[0] = vec4(tex, 0, 1);

//  vec4 tang = vec4(gl_NormalMatrix * tangent.xyz, 0);//vec4(tangent, 0); // tagent in eye sapce
//  vec4 tang = vec4(gl_NormalMatrix * gl_Color.xyz, 0);//vec4(tangent, 0); // tagent in eye sapce
  vec3 tang = gl_NormalMatrix * vec3(0, 0, 1);//vec4(tangent, 0); // tagent in eye sapce
  tang = tang - dot(tang, normalEye) * normalEye;
//  tang[3] = 0.0;
  tang = normalize(tang);
  vec3 binormal = normalize(cross(normalEye.xyz, tang.xyz));
  // eye space to tagent space transform
  mat3 toObjectLocal = mat3(tang.x, binormal.x, normalEye.x,
                            tang.y, binormal.y, normalEye.y,
                            tang.z, binormal.z, normalEye.z);
  LightDir.xyz = normalize(toObjectLocal * (-positionEye.xyz)); //vec3(0,0,0)
  ViewDir = LightDir;
}
